"""ICE database: plasma fraction unbound used by the ICE IVIVE tool (NICEATM Integrated Chemical Environment).

Assays: httk Human and Rat Plasma Fraction Unbound (measured values in EPA's httk package) and
OPERA Human Plasma Fraction Unbound (QSAR prediction with its applicability domain).
Compounds are sent by InChIKey (25 per request), then any compound still without a record is sent alone
by SMILES. ICE answers a chemid it cannot read with its whole database, so a reply over 2 MB is treated
as "not found". A record is kept only when its InChIKey (qsarReadyId) equals the compound's, or shares
the first block (same skeleton, stereochemistry may differ).

Input : table compounds; ICE REST API (replies cached in output/cache/ice)
Output: table ice_fraction_unbound, and output/ice.xlsx
Usage : python Source_Code/ice.py <database.db>
"""
import hashlib
import json
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

API = "https://ice.ntp.niehs.nih.gov/api/v1/search"
ASSAYS = ["httk, Human Plasma Fraction Unbound", "httk, Rat Plasma Fraction Unbound",
          "OPERA, Human Plasma Fraction Unbound"]
LIMIT = 2_000_000   # bytes; a real reply is a few kB
db = Path(sys.argv[1])
CACHE = db.parent / "cache" / "ice"


def post(chemids):
    """POST once; the reply is saved in the cache and re-used on every later run."""
    body = {"chemids": chemids, "assays": ASSAYS}
    path = CACHE / (hashlib.sha1(json.dumps(body, sort_keys=True).encode()).hexdigest() + ".json")
    if path.exists():
        return json.loads(path.read_text())
    for wait in (10, 20, 40, 80, 120, 120):
        try:
            r = requests.post(API, json=body, timeout=180, stream=True)
            if r.status_code == 200:
                data = b""
                for chunk in r.iter_content(1 << 16):
                    data += chunk
                    if len(data) > LIMIT:
                        break
                r.close()
                records = json.loads(data).get("endPoints", []) if len(data) <= LIMIT else []
                break
            print(f"    ICE HTTP {r.status_code}, retrying in {wait} s")
        except requests.RequestException as e:
            print(f"    ICE {type(e).__name__}, retrying in {wait} s")
        time.sleep(wait)
    else:
        sys.exit("ICE request failed")
    saved = {"retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "records": records}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(saved))
    time.sleep(2)   # ICE limits the request rate
    return saved


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


con = sqlite3.connect(db)
compounds = pd.read_sql("SELECT compound_id, name, smiles_input, inchikey FROM compounds "
                        "WHERE inchikey IS NOT NULL ORDER BY compound_id", con)
exact = dict(zip(compounds.inchikey, compounds.compound_id))
skeleton = compounds.groupby(compounds.inchikey.str[:14]).compound_id.apply(list).to_dict()


def matched(replies):
    """Database rows for the ICE records that belong to one of the compounds."""
    rows = set()
    for rec in (rec for reply in replies for rec in reply["records"]):
        key = rec.get("qsarReadyId") or ""
        if rec.get("assay") not in ASSAYS or len(key) < 14:
            continue
        ids, method = (([exact[key]], "inchikey") if key in exact
                       else (skeleton.get(key[:14], []), "inchikey_connectivity"))
        rows |= {(cid, rec["assay"], rec.get("endpoint") or "", num(rec.get("value")), rec.get("unit") or None,
                  rec.get("species") or ("Rat" if "Rat" in rec["assay"] else "Human"), rec.get("reference") or None,
                  rec.get("substanceName") or None, rec.get("dtxsid") or None, rec.get("casrn") or None, key, method)
                 for cid in ids}
    return rows


replies = [post(list(compounds.inchikey[i:i + 25])) for i in range(0, len(compounds), 25)]
missing = compounds[~compounds.compound_id.isin({r[0] for r in matched(replies)})]
print(f"    ICE: {len(compounds) - len(missing)} found by InChIKey, trying {len(missing)} by SMILES")
replies += [post([smiles]) for smiles in missing.smiles_input]
rows = matched(replies)
con.executemany("INSERT INTO ice_fraction_unbound (compound_id, assay, endpoint, value, unit, species, reference, "
                "substance_name, dtxsid, casrn, ice_inchikey, match_method) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                sorted(rows, key=lambda r: tuple("" if v is None else str(v) for v in r)))
retrieved = min(reply["retrieved_utc"] for reply in replies)
con.execute("INSERT INTO sources VALUES (?, ?, ?, ?)", ("ICE", "httk / OPERA versions in 'reference'", retrieved, API))
con.commit()

records = pd.read_sql("SELECT c.name, f.* FROM ice_fraction_unbound f JOIN compounds c USING (compound_id) "
                      "ORDER BY f.compound_id, f.assay, f.endpoint", con)
# Summary: exact InChIKey records when a compound has them, otherwise the same-skeleton records.
has_exact = records.groupby(["compound_id", "assay"]).match_method.transform(lambda m: (m == "inchikey").any())
best = records[~has_exact | (records.match_method == "inchikey")]
summary = compounds[["compound_id", "name"]]
for assay, label in zip(ASSAYS, ["httk_human", "httk_rat", "opera_human"]):
    part = best[(best.assay == assay) & (best.endpoint == "Fu")].groupby("compound_id").agg(
        **{f"{label}_fu": ("value", lambda v: "; ".join(map(str, v))),
           f"{label}_reference": ("reference", lambda v: "; ".join(map(str, v.dropna()))),
           f"{label}_match": ("match_method", lambda v: "; ".join(sorted(set(v))))})
    summary = summary.merge(part, on="compound_id", how="left")
domain = best[(best.assay == ASSAYS[2]) & (best.endpoint == "Applicability_Domain")]
summary = summary.merge(domain.groupby("compound_id").value.max().rename("opera_in_domain (1 = yes)"),
                        on="compound_id", how="left").drop(columns="opera_human_reference")

with pd.ExcelWriter(db.parent / "ice.xlsx") as xl:
    summary.to_excel(xl, sheet_name="Fraction unbound", index=False)
    records.to_excel(xl, sheet_name="All records", index=False)
    compounds[~compounds.compound_id.isin(records.compound_id)].to_excel(xl, sheet_name="Not found", index=False)
    pd.DataFrame([["API", API], ["Assays", "; ".join(ASSAYS)], ["First downloaded (UTC)", retrieved],
                  ["Searched by", "InChIKey (25 per request), then SMILES one at a time for compounds still missing"],
                  ["match_method", "inchikey = exact; inchikey_connectivity = same skeleton, stereo may differ"]],
                 columns=["Item", "Value"]).to_excel(xl, sheet_name="Source", index=False)
print(f"ice: {records.compound_id.nunique()}/{len(compounds)} compounds with fraction unbound -> ice.xlsx")
