"""ChEMBL database: biological targets and individual ADME parameters (human), including protein binding.

Each compound is matched to one ChEMBL molecule: exact InChIKey, then same skeleton (first InChIKey
block), then name, then synonym; salt forms are replaced by their parent. For that molecule:
  mechanisms of action + targets                 -> chembl_mechanisms
  each ADME parameter in ADME, human assays only -> chembl_adme (sheets per parameter and in vivo / in vitro)
  binding records with pChEMBL (assay type B)    -> chembl_binding
Errors (SD, SEM, CI) are copied only where ChEMBL reports them; nothing is calculated.
Each ADME value is labelled from ChEMBL's assay format: organism-based = in vivo, no format given = unspecified,
any other format (microsome, cell, cell-free, tissue, protein) = in vitro.

Input : table compounds; ChEMBL web services (replies cached in output/cache/chembl)
Output: tables chembl_molecules, chembl_mechanisms, chembl_adme, chembl_binding, and output/chembl.xlsx
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

ADME = {"Absorption": ["F", "Fa", "Papp", "Tmax", "Cmax", "AUC"],
        "Distribution": ["Fu", "PPB", "Vdss", "Vd"],
        "Metabolism": ["CL", "CLH", "T1/2", "MRT", "Stability", "Km", "Vmax"],
        "Excretion": ["CL_renal", "Drug excretion"]}
GROUP = {parameter: group for group, parameters in ADME.items() for parameter in parameters}
SYSTEMS = ["in vivo", "in vitro", "unspecified"]
ERRORS = {"SD", "SEM", "SE", "S.D.", "S.E.M.", "CI", "CI_LOWER", "CI_UPPER", "CV", "ERROR"}
API = "https://www.ebi.ac.uk/chembl/api/data"
MOLECULE = "molecule_chembl_id,pref_name,max_phase,molecule_hierarchy"
ADME_FIELDS = ("activity_id,assay_chembl_id,assay_description,standard_type,standard_relation,standard_value,"
               "standard_upper_value,standard_units,document_chembl_id,activity_properties")
BINDING_FIELDS = ("activity_id,assay_chembl_id,assay_description,standard_type,standard_relation,standard_value,"
                  "standard_units,pchembl_value,target_chembl_id,target_pref_name,target_organism,document_chembl_id")
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
        return c.compound_id, m, method, [], [], []
    return (c.compound_id, m, method,
            get_all("mechanism", "mechanisms", molecule_chembl_id=cid),
            get_all("activity", "activities", molecule_chembl_id=cid, assay_type="A",
                    standard_type__in=",".join(GROUP), only=ADME_FIELDS),
            get_all("activity", "activities", molecule_chembl_id=cid, assay_type="B",
                    pchembl_value__isnull="false", only=BINDING_FIELDS))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def errors(record):
    """Error terms exactly as ChEMBL reports them in the record's properties."""
    found = []
    for p in record.get("activity_properties") or []:
        kind = (p.get("standard_type") or p.get("type") or "").upper()
        if kind in ERRORS or "DEVIATION" in kind or "ERROR" in kind or "CONFIDENCE" in kind:
            value = p.get("standard_value") or p.get("value") or p.get("standard_text_value") or ""
            found.append(" ".join(str(s) for s in (kind, p.get("standard_relation") or "", value,
                                                   p.get("standard_units") or p.get("units") or "") if s))
    return "; ".join(found) or None


# Study details: column -> ChEMBL property / assay parameter names, copied as reported.
DETAILS = {"dosed_drug": ["DOSED_COMPOUND_NAME"], "interacting_drug": ["INTERACTING_COMPOUND"],
           "health_status": ["HEALTH_STATUS"], "gender": ["GENDER"], "n_subjects": ["N"],
           "age_range": ["AGE_RANGE"], "dietary_status": ["DIETARY_STATUS"], "dose": ["DOSE"], "route": ["ROUTE"],
           "regimen": ["REGIMEN"], "formulation": ["FORMULATION", "DOSAGE_FORM"], "sample_tissue": ["TISSUE"]}


def details(record, assay):
    """Study details from the value's properties, else from its assay's parameters."""
    def read(items):
        out = {}
        for p in items or []:
            kind = (p.get("standard_type") or p.get("type") or "").upper()
            text = p.get("standard_text_value") or p.get("text_value")
            if text is None:
                value = num(p.get("standard_value") if p.get("standard_value") is not None else p.get("value"))
                text = None if value is None else f"{value:g} {p.get('standard_units') or p.get('units') or ''}".strip()
            if text:
                out.setdefault(kind, []).append(str(text))
        return out
    own, study = read(record.get("activity_properties")), read(assay.get("assay_parameters"))
    row = {}
    for column, kinds in DETAILS.items():
        values = next((own[k] for k in kinds if k in own), None) or next((study[k] for k in kinds if k in study), [])
        row[column] = "; ".join(dict.fromkeys(values)) or None
    row["sample_tissue"] = row["sample_tissue"] or assay.get("assay_tissue")
    return row


con = sqlite3.connect(db)
compounds = pd.read_sql("SELECT compound_id, name, inchikey FROM compounds ORDER BY compound_id", con)
status = get("status")
results = []
with ThreadPoolExecutor(max_workers=4) as pool:
    for result in pool.map(fetch, compounds.itertuples(index=False)):
        results.append(result)
        if len(results) % 25 == 0:
            print(f"    ChEMBL: {len(results)}/{len(compounds)} compounds")


def lookup(endpoint, key, id_field, ids, fields):
    ids = sorted(ids)
    return {x[id_field]: x for i in range(0, len(ids), 50)
            for x in get_all(endpoint, key, only=fields, **{f"{id_field}__in": ",".join(ids[i:i + 50])})}


targets = lookup("target", "targets", "target_chembl_id",
                 {x["target_chembl_id"] for r in results for x in r[3] if x.get("target_chembl_id")},
                 "target_chembl_id,pref_name,target_type,organism")
assays = lookup("assay", "assays", "assay_chembl_id", {a["assay_chembl_id"] for r in results for a in r[4]},
                "assay_chembl_id,assay_organism,assay_tax_id,bao_label,assay_tissue,assay_parameters")


def system(fmt):
    return "in vivo" if fmt == "organism-based format" else "unspecified" if fmt in (None, "assay format") \
        else "in vitro"

for compound_id, m, method, mechanisms, adme, binding in results:
    cid = m.get("molecule_chembl_id")
    con.execute("INSERT INTO chembl_molecules VALUES (?, ?, ?, ?, ?)",
                (compound_id, cid, m.get("pref_name"), num(m.get("max_phase")), method))
    for x in mechanisms:
        t = targets.get(x.get("target_chembl_id"), {})
        con.execute("INSERT INTO chembl_mechanisms VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (compound_id, x["mec_id"], cid, x.get("mechanism_of_action"), x.get("action_type"),
                     x.get("target_chembl_id"), t.get("pref_name"), t.get("target_type"), t.get("organism")))
    for a in {a["activity_id"]: a for a in adme}.values():
        assay = assays.get(a["assay_chembl_id"], {})
        if str(assay.get("assay_tax_id")) != "9606" or a["standard_type"] not in GROUP:   # human assays only
            continue
        row = {"compound_id": compound_id, "activity_id": a["activity_id"], "molecule_chembl_id": cid,
               "adme_group": GROUP[a["standard_type"]], "parameter": a["standard_type"],
               "relation": a.get("standard_relation"), "value": num(a.get("standard_value")),
               "upper_value": num(a.get("standard_upper_value")), "units": a.get("standard_units"),
               "error_reported": errors(a), "assay_chembl_id": a.get("assay_chembl_id"),
               "assay_description": a.get("assay_description"), "assay_organism": assay.get("assay_organism"),
               "assay_format": assay.get("bao_label"), "system": system(assay.get("bao_label")),
               **details(a, assay), "document_chembl_id": a.get("document_chembl_id")}
        con.execute(f"INSERT INTO chembl_adme ({', '.join(row)}) VALUES ({', '.join('?' * len(row))})",
                    list(row.values()))
    for a in {a["activity_id"]: a for a in binding}.values():
        con.execute("INSERT INTO chembl_binding VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (compound_id, a["activity_id"], cid, a.get("assay_chembl_id"), a.get("assay_description"),
                     a.get("standard_type"), a.get("standard_relation"), num(a.get("standard_value")),
                     a.get("standard_units"), num(a.get("pchembl_value")), a.get("target_chembl_id"),
                     a.get("target_pref_name"), a.get("target_organism"), a.get("document_chembl_id")))
release = f"{status['chembl_db_version']} ({status['chembl_release_date']})"
retrieved = json.loads(cached("status", {}).read_text())["retrieved_utc"]
con.execute("INSERT INTO sources VALUES (?, ?, ?, ?)", ("ChEMBL", release, retrieved, API))
con.commit()

# Spreadsheet: summary (every value listed), one sheet per ADME parameter, targets, binding.
adme = pd.read_sql("SELECT c.name, a.* FROM chembl_adme a JOIN compounds c USING (compound_id) "
                   "ORDER BY a.compound_id, a.activity_id", con)


def shown(r):
    relation = "" if r.relation in (None, "=") else r.relation
    upper = f"-{r.upper_value:g}" if pd.notna(r.upper_value) else ""
    value = f"{r.value:g}" if pd.notna(r.value) else "?"
    return f"{relation}{value}{upper} {r.units or ''}".strip()


adme["shown"] = [shown(r) for r in adme.itertuples(index=False)]


def summary(rows):
    """Compounds x parameters, every value listed."""
    table = rows.groupby(["compound_id", "parameter"])["shown"].agg("; ".join).unstack()
    table = compounds[["compound_id", "name"]].merge(
        table.reindex(columns=list(GROUP)).reset_index(), on="compound_id", how="left")
    table.columns = ["compound_id", "name"] + [f"{GROUP[p]}: {p}" for p in GROUP]
    return table


named = "SELECT c.name, t.* FROM {} t JOIN compounds c USING (compound_id) ORDER BY t.compound_id"
counts = ", ".join(f"(SELECT COUNT(*) FROM chembl_adme x WHERE x.compound_id = c.compound_id AND x.system = '{s}') "
                   f"AS \"n_adme_{s}\"" for s in SYSTEMS)
sheets = {"Molecules": pd.read_sql(
    "SELECT c.compound_id, c.name, c.inchikey, m.molecule_chembl_id, m.pref_name, m.max_phase, m.match_method, "
    "(SELECT COUNT(*) FROM chembl_mechanisms x WHERE x.compound_id = c.compound_id) AS n_mechanisms, "
    f"{counts}, (SELECT COUNT(*) FROM chembl_binding x WHERE x.compound_id = c.compound_id) AS n_binding "
    "FROM compounds c JOIN chembl_molecules m USING (compound_id) ORDER BY c.compound_id", con)}
for s in SYSTEMS:
    sheets[f"ADME {s}"] = summary(adme[adme.system == s])
for p in GROUP:
    for s in SYSTEMS:
        rows = adme[(adme.parameter == p) & (adme.system == s)]
        if len(rows):
            sheets[f"{GROUP[p][0]} {p} {s}".replace("/", "-")] = rows.drop(columns="shown")
sheets["Mechanism targets"] = pd.read_sql(named.format("chembl_mechanisms"), con)
sheets["Binding target summary"] = pd.read_sql(named.format("chembl_target_summary"), con)
sheets["Binding (pChEMBL)"] = pd.read_sql(named.format("chembl_binding"), con)
sheets["Source"] = pd.DataFrame(
    [["ChEMBL release", release], ["First downloaded (UTC)", retrieved], ["Web services", API],
     ["Matching", "exact InChIKey > same skeleton > name > synonym; '_parent' = salt replaced by parent"],
     ["ADME parameters", "; ".join(f"{g}: {', '.join(p)}" for g, p in ADME.items())],
     ["ADME species", "human only (assay taxonomy 9606)"],
     ["In vivo / in vitro", "from ChEMBL assay format: organism-based = in vivo; microsome, cell, cell-free, tissue, "
                            "protein = in vitro; not given = unspecified (column assay_format keeps the original)"],
     ["Errors", "SD / SEM / CI copied only where ChEMBL reports them; none calculated"],
     ["Study details", "dosed_drug, interacting_drug, health_status, gender, n_subjects, age_range, dietary_status, "
                       "dose, route, regimen, formulation, sample_tissue: as reported by ChEMBL for the value, "
                       "else for its assay; blank when not reported"],
     ["Sheet names", "A = Absorption, D = Distribution, M = Metabolism, E = Excretion"],
     ["Binding", "assay type B records that have a pChEMBL value"]], columns=["Item", "Value"])
with pd.ExcelWriter(db.parent / "chembl.xlsx") as xl:
    for sheet, df in sheets.items():
        df.to_excel(xl, sheet_name=sheet, index=False)
matched = sum(1 for r in results if r[1])
split = ", ".join(f"{(adme.system == s).sum()} {s}" for s in SYSTEMS)
print(f"chembl: {release}, {matched}/{len(results)} compounds matched, human ADME values: {split} -> chembl.xlsx")
