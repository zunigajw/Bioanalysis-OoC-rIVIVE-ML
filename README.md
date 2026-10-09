# Bioanalysis-OoC-rIVIVE-ML

Database Curation: one command builds a database of physicochemical, ADME and protein-binding data for the
compounds in `Project_Files/compound_list.csv`, plus one spreadsheet per data source.

```
git clone https://github.com/zunigajw/Bioanalysis-OoC-rIVIVE-ML && cd Bioanalysis-OoC-rIVIVE-ML && just
```

Requires `git`, `podman` and `just` (see the course setup). Everything else runs inside the container.

## Pipeline flow

```mermaid
flowchart TD
    J["<b>justfile</b> (Database Curation)<br/>the one command: just"]
    C["<b>Containerfile</b> (Podman)<br/>builds the same Linux environment on any computer:<br/>Python 3.12, RDKit, pandas, pinned versions"]
    G["<b>GenDB</b> (makefile)<br/>runs every step in order"]
    S["<b>Data_Pipeline.sql</b><br/>creates the empty database tables"]
    P[/"<b>Project_Files/compound_list.csv</b><br/>128 compounds: name + SMILES"/]

    R["<b>rdkit_features.py</b> (RDKit)<br/>canonical SMILES and InChIKey<br/>217 descriptors: MolWt, LogP, TPSA,<br/>H-bond donors/acceptors ...<br/>Morgan fingerprints: radius 2, 2048 bits"]
    M["<b>chembl.py</b> (ChEMBL 37)<br/>mechanism-of-action targets<br/>human ADME, split in vivo / in vitro / unspecified:<br/>Absorption: F, Fa, Papp, Tmax, Cmax, AUC<br/>Distribution: Fu, PPB, Vdss, Vd<br/>Metabolism: CL, CLH, T1/2, MRT, stability, Km, Vmax<br/>Excretion: CL_renal, drug excretion<br/>study details: dosed drug, health status, gender,<br/>N, age, dose, route, formulation<br/>binding affinities (pChEMBL)"]
    I["<b>ice.py</b> (NICEATM ICE)<br/>plasma fraction unbound:<br/>httk human and rat (measured),<br/>OPERA human (predicted)"]

    XR[/"<b>rdkit.xlsx</b><br/>database tables: compounds,<br/>rdkit_descriptors, morgan_fingerprints"/]
    XM[/"<b>chembl.xlsx</b><br/>database tables: chembl_molecules,<br/>chembl_mechanisms, chembl_adme, chembl_binding"/]
    XI[/"<b>ice.xlsx</b><br/>database table: ice_fraction_unbound"/]
    D[/"<b>ooc_drug_ml.db</b><br/>one SQLite file holding all of the tables above,<br/>linked by compound_id"/]

    J --> C --> G --> S --> P
    P --> R & M & I
    R --> XR
    M --> XM
    I --> XI
    XR & XM & XI --> D
```

Rectangles are code; slanted boxes are data files (the input list and every output).

## What `just` creates in `output/` (not stored in the repository)

| File | What is in it |
|---|---|
| `ooc_drug_ml.db` | The database. Tables: `compounds`, `rdkit_descriptors`, `morgan_fingerprints`, `chembl_molecules`, `chembl_mechanisms`, `chembl_adme`, `chembl_binding`, `ice_fraction_unbound`, `sources`; columns are described in `Data_Pipeline.sql` |
| `rdkit.xlsx` | Compounds (SMILES, InChIKey), the 217 RDKit descriptors, Morgan fingerprint bits |
| `chembl.xlsx` | ChEMBL match per compound; ADME summaries (in vivo, in vitro, unspecified); one sheet per ADME parameter and system with every human value, its units, reported error and study details; mechanism targets; binding |
| `ice.xlsx` | Plasma fraction unbound per compound (httk human and rat, OPERA predicted with applicability domain) |
