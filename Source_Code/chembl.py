"""ChEMBL database: biological targets, ADME properties and protein-binding parameters.

Each compound is matched to one ChEMBL molecule: exact InChIKey, then same skeleton (first InChIKey
block), then name, then synonym; salt forms are replaced by their parent. For that molecule:
  mechanisms of action + targets          -> chembl_mechanisms
  ADME records (assay type A)             -> chembl_activities, assay_type 'A'
  binding records with pChEMBL (type B)   -> chembl_activities, assay_type 'B'
Protein-binding parameters (Fu, PPB) are the ADME records in view chembl_protein_binding.

Input : table compounds; ChEMBL web services (replies cached in output/cache/chembl)
Output: tables chembl_molecules, chembl_mechanisms, chembl_activities, and output/chembl.xlsx
Usage : python Source_Code/chembl.py <database.db>
"""
import hashlib
import json
import sqlite3
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

API = "https://www.ebi.ac.uk/chembl/api/data"
MOLECULE = "molecule_chembl_id,pref_name,max_phase,molecule_hierarchy"
ACTIVITY = ("activity_id,assay_chembl_id,assay_description,assay_type,standard_type,standard_relation,"
            "standard_value,standard_units,pchembl_value,target_chembl_id,target_pref_name,"
            "target_organism,document_chembl_id")
db = Path(sys.argv[1])
CACHE = db.parent / "cache" / "chembl"


def cached(endpoint, params):
    return CACHE / (hashlib.sha1(json.dumps([endpoint, params], sort_keys=True).encode()).hexdigest() + ".json")


def get(endpoint, **params):
    """GET once; the reply is saved in the cache and re-used on every later run."""
    path = cached(endpoint, params)
    if path.exists():
        return json.loads(path.read_text())["reply"]
    for wait in (5, 10, 20, 40, 80, 120):
        try:
            r = requests.get(f"{API}/{endpoint}.json", params=params, timeout=180)
            if r.status_code == 200:
                break
            print(f"    ChEMBL HTTP {r.status_code}, retrying in {wait} s")
        except requests.RequestException as e:
            print(f"    ChEMBL {type(e).__name__}, retrying in {wait} s")
        time.sleep(wait)
    else:
        sys.exit(f"ChEMBL request failed: {endpoint} {params}")
    path.parent.mkdir(parents=True, exist_ok=True)
    utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    path.write_text(json.dumps({"retrieved_utc": utc, "reply": r.json()}))
    return r.json()


def get_all(endpoint, key, **params):
    """All records of a list query, following ChEMBL's pages of 1000."""
    records, offset = [], 0
    while True:
        page = get(endpoint, **params, limit=1000, offset=offset)
        records += page[key]
        if not page["page_meta"]["next"]:
            return records
        offset += 1000


def parent(m):
    return (m.get("molecule_hierarchy") or {}).get("parent_chembl_id") or m["molecule_chembl_id"]


def match(c):
    searches = [("name", {"pref_name__iexact": c.name}),
                ("synonym", {"molecule_synonyms__molecule_synonym__iexact": c.name})]
    if c.inchikey:
        searches = [("inchikey", {"molecule_structures__standard_inchi_key": c.inchikey}),
                    ("inchikey_connectivity",
                     {"molecule_structures__standard_inchi_key__startswith": c.inchikey[:14]})] + searches
    for method, query in searches:
        found = get_all("molecule", "molecules", only=MOLECULE, **query)
        if found:
            # prefer a parent molecule, then the highest clinical phase, then the oldest ChEMBL ID
            m = min(found, key=lambda m: (parent(m) != m["molecule_chembl_id"], -float(m.get("max_phase") or -1),
                                          int(m["molecule_chembl_id"][6:])))
            if parent(m) != m["molecule_chembl_id"]:
                m, method = get_all("molecule", "molecules", only=MOLECULE, molecule_chembl_id=parent(m))[0], \
                    method + "_parent"
            return m, method
    return {}, "not_found"


def fetch(c):
    m, method = match(c)
    cid = m.get("molecule_chembl_id")
    if not cid:
        return c.compound_id, m, method, [], []
    mechanisms = get_all("mechanism", "mechanisms", molecule_chembl_id=cid)
    activities = (get_all("activity", "activities", molecule_chembl_id=cid, assay_type="A", only=ACTIVITY)
                  + get_all("activity", "activities", molecule_chembl_id=cid, assay_type="B",
                            pchembl_value__isnull="false", only=ACTIVITY))
    return c.compound_id, m, method, mechanisms, activities


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


con = sqlite3.connect(db)
compounds = pd.read_sql("SELECT compound_id, name, inchikey FROM compounds ORDER BY compound_id", con)
status = get("status")
results = []
with ThreadPoolExecutor(max_workers=4) as pool:
    for result in pool.map(fetch, compounds.itertuples(index=False)):
        results.append(result)
        if len(results) % 25 == 0:
            print(f"    ChEMBL: {len(results)}/{len(compounds)} compounds")

target_ids = sorted({x["target_chembl_id"] for r in results for x in r[3] if x.get("target_chembl_id")})
targets = {t["target_chembl_id"]: t for i in range(0, len(target_ids), 50)
           for t in get_all("target", "targets", target_chembl_id__in=",".join(target_ids[i:i + 50]),
                            only="target_chembl_id,pref_name,target_type,organism")}

for compound_id, m, method, mechanisms, activities in results:
    cid = m.get("molecule_chembl_id")
    con.execute("INSERT INTO chembl_molecules VALUES (?, ?, ?, ?, ?)",
                (compound_id, cid, m.get("pref_name"), num(m.get("max_phase")), method))
    for x in mechanisms:
        t = targets.get(x.get("target_chembl_id"), {})
        con.execute("INSERT INTO chembl_mechanisms VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (compound_id, x["mec_id"], cid, x.get("mechanism_of_action"), x.get("action_type"),
                     x.get("target_chembl_id"), t.get("pref_name"), t.get("target_type"), t.get("organism")))
    for a in {a["activity_id"]: a for a in activities}.values():
        con.execute("INSERT INTO chembl_activities VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (compound_id, a["activity_id"], cid, a["assay_type"], a.get("assay_chembl_id"),
                     a.get("assay_description"), a.get("standard_type"), a.get("standard_relation"),
                     num(a.get("standard_value")), a.get("standard_units"), num(a.get("pchembl_value")),
                     a.get("target_chembl_id"), a.get("target_pref_name"), a.get("target_organism"),
                     a.get("document_chembl_id")))
release = f"{status['chembl_db_version']} ({status['chembl_release_date']})"
retrieved = json.loads(cached("status", {}).read_text())["retrieved_utc"]
con.execute("INSERT INTO sources VALUES (?, ?, ?, ?)", ("ChEMBL", release, retrieved, API))
con.commit()

named = "SELECT c.name, t.* FROM ({}) t JOIN compounds c USING (compound_id) ORDER BY t.compound_id"
counts = ", ".join(f"(SELECT COUNT(*) FROM {v} x WHERE x.compound_id = c.compound_id) AS n_{v[7:]}"
                   for v in ("chembl_mechanisms", "chembl_adme", "chembl_protein_binding", "chembl_binding"))
sheets = {
    "Molecules": f"SELECT c.compound_id, c.name, c.inchikey, m.molecule_chembl_id, m.pref_name, m.max_phase, "
                 f"m.match_method, {counts} FROM compounds c JOIN chembl_molecules m USING (compound_id) "
                 f"ORDER BY c.compound_id",
    "Mechanism targets": named.format("SELECT * FROM chembl_mechanisms"),
    "ADME": named.format("SELECT * FROM chembl_adme"),
    "Protein binding": named.format("SELECT * FROM chembl_protein_binding"),
    "Binding target summary": named.format("SELECT * FROM chembl_target_summary"),
    "Binding (pChEMBL)": named.format("SELECT * FROM chembl_binding"),
}
with pd.ExcelWriter(db.parent / "chembl.xlsx") as xl:
    for sheet, sql in sheets.items():
        pd.read_sql(sql, con).to_excel(xl, sheet_name=sheet, index=False)
    pd.DataFrame([["ChEMBL release", release], ["First downloaded (UTC)", retrieved], ["Web services", API],
                  ["Matching", "exact InChIKey > same skeleton > name > synonym; '_parent' = salt replaced by parent"],
                  ["Protein binding", "ADME records with standard_type Fu or PPB"],
                  ["Binding", "assay type B records that have a pChEMBL value"]],
                 columns=["Item", "Value"]).to_excel(xl, sheet_name="Source", index=False)
matched = sum(1 for r in results if r[1])
print(f"chembl: {release}, {matched}/{len(results)} compounds matched -> chembl.xlsx")
