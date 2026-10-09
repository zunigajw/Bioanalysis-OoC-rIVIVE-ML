"""RDKit database: compound list -> compounds, 2D descriptors and Morgan fingerprints.

Input : Project_Files/compound_list.csv
Output: tables compounds, rdkit_descriptors, morgan_fingerprints, and output/rdkit.xlsx
Usage : python Source_Code/rdkit_features.py <database.db>
"""
import hashlib
import math
import sqlite3
import sys
from pathlib import Path

import pandas as pd
import rdkit
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
SRC = Path(__file__).resolve().parent.parent / "Project_Files" / "compound_list.csv"
RADIUS, N_BITS = 2, 2048   # same Morgan settings as qssr\Morgan.csv
NAMES = [name for name, _ in Descriptors.descList]


def number(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


db = Path(sys.argv[1])
con = sqlite3.connect(db)
morgan = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=N_BITS)

compounds, descriptors, wide, fingerprints, bits, first_id = [], [], [], [], [], {}
for cid, r in enumerate(pd.read_csv(SRC).itertuples(index=False), 1):
    name, smiles = str(r.Compound).strip(), str(r.SMILES).strip()
    mol = Chem.MolFromSmiles(smiles)
    key = Chem.MolToInchiKey(mol) if mol else None
    compounds.append((cid, name, smiles, Chem.MolToSmiles(mol) if mol else None, key,
                      int(mol is not None), first_id.get(key) if key else None))
    if key:
        first_id.setdefault(key, cid)
    if mol is None:
        continue
    values = {n: number(v) for n, v in Descriptors.CalcMolDescriptors(mol).items()}
    descriptors += [(cid, n, values.get(n)) for n in NAMES]
    wide.append({"compound_id": cid, "name": name, **{n: values.get(n) for n in NAMES}})
    fp = morgan.GetFingerprint(mol)
    fingerprints.append((cid, RADIUS, N_BITS, " ".join(map(str, fp.GetOnBits())), fp.ToBitString()))
    bits.append([cid, name] + [int(b) for b in fp.ToBitString()])

con.executemany("INSERT INTO compounds VALUES (?, ?, ?, ?, ?, ?, ?)", compounds)
con.executemany("INSERT INTO rdkit_descriptors VALUES (?, ?, ?)", descriptors)
con.executemany("INSERT INTO morgan_fingerprints VALUES (?, ?, ?, ?, ?)", fingerprints)
checksum = "sha256:" + hashlib.sha256(SRC.read_bytes()).hexdigest()[:16]
con.executemany("INSERT INTO sources VALUES (?, ?, ?, ?)", [
    ("compound_list", checksum, None, "Project_Files/compound_list.csv"),
    ("RDKit", rdkit.__version__, None, "Descriptors.descList; Morgan radius 2, 2048 bits")])
con.commit()

source = pd.DataFrame([["Input", "Project_Files/compound_list.csv"], ["Input checksum", checksum],
                       ["RDKit version", rdkit.__version__], ["Descriptors", len(NAMES)],
                       ["Morgan radius / bits", f"{RADIUS} / {N_BITS}"],
                       ["Compounds", len(compounds)], ["SMILES not parsed", len(compounds) - len(wide)]],
                      columns=["Item", "Value"])
sheets = {"Compounds": pd.read_sql("SELECT * FROM compounds ORDER BY compound_id", con),
          "Descriptors": pd.DataFrame(wide),
          "Morgan": pd.DataFrame(bits, columns=["compound_id", "name"] + list(map(str, range(N_BITS)))),
          "Source": source}
with pd.ExcelWriter(db.parent / "rdkit.xlsx") as xl:
    for sheet, df in sheets.items():
        df.to_excel(xl, sheet_name=sheet, index=False)
print(f"rdkit_features: {len(compounds)} compounds, {len(NAMES)} descriptors, Morgan {N_BITS} bits -> rdkit.xlsx")
