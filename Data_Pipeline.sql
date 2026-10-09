-- Data Pipeline: database schema (SQLite), created empty by GenDB before the scripts run.
PRAGMA foreign_keys = ON;

-- Source and version of each block of data.
CREATE TABLE sources (
  source        TEXT PRIMARY KEY,   -- compound_list, RDKit, ChEMBL, ICE
  version       TEXT,
  retrieved_utc TEXT,
  location      TEXT
);

-- rdkit_features.py: Project_Files/compound_list.csv standardised with RDKit.
CREATE TABLE compounds (
  compound_id      INTEGER PRIMARY KEY,                         -- row number in the compound list
  name             TEXT NOT NULL,
  smiles_input     TEXT NOT NULL,
  canonical_smiles TEXT,
  inchikey         TEXT,                                        -- used to search ChEMBL and ICE
  parse_ok         INTEGER NOT NULL CHECK (parse_ok IN (0, 1)),
  duplicate_of     INTEGER REFERENCES compounds (compound_id)   -- earlier ID with the same InChIKey
);

CREATE TABLE rdkit_descriptors (
  compound_id INTEGER NOT NULL REFERENCES compounds (compound_id),
  descriptor  TEXT NOT NULL,
  value       REAL,
  PRIMARY KEY (compound_id, descriptor)
);

CREATE TABLE morgan_fingerprints (
  compound_id INTEGER PRIMARY KEY REFERENCES compounds (compound_id),
  radius      INTEGER NOT NULL,
  n_bits      INTEGER NOT NULL,
  on_bits     TEXT NOT NULL,   -- space-separated indices of the bits set to 1
  bitstring   TEXT NOT NULL    -- n_bits characters of 0/1, bit 0 first
);

-- chembl.py: matched molecule, mechanism targets, individual ADME parameters and binding.
CREATE TABLE chembl_molecules (
  compound_id        INTEGER PRIMARY KEY REFERENCES compounds (compound_id),
  molecule_chembl_id TEXT,
  pref_name          TEXT,
  max_phase          REAL,
  match_method       TEXT NOT NULL   -- inchikey, inchikey_connectivity, name, synonym, not_found
);

CREATE TABLE chembl_mechanisms (
  compound_id         INTEGER NOT NULL REFERENCES compounds (compound_id),
  mec_id              INTEGER NOT NULL,
  molecule_chembl_id  TEXT NOT NULL,
  mechanism_of_action TEXT,
  action_type         TEXT,
  target_chembl_id    TEXT,
  target_pref_name    TEXT,
  target_type         TEXT,
  target_organism     TEXT,
  PRIMARY KEY (compound_id, mec_id)
);

-- Individual ADME parameters (human assays only), one row per reported value.
CREATE TABLE chembl_adme (
  compound_id        INTEGER NOT NULL REFERENCES compounds (compound_id),
  activity_id        INTEGER NOT NULL,
  molecule_chembl_id TEXT NOT NULL,
  adme_group         TEXT NOT NULL CHECK (adme_group IN ('Absorption', 'Distribution', 'Metabolism', 'Excretion')),
  parameter          TEXT NOT NULL,   -- ChEMBL standard_type, e.g. F, Fu, CL, T1/2, CL_renal
  relation           TEXT,            -- =, <, >, ~
  value              REAL,
  upper_value        REAL,            -- upper end when a range was reported
  units              TEXT,
  error_reported     TEXT,            -- SD / SEM / CI exactly as reported in ChEMBL; never calculated
  assay_chembl_id    TEXT,
  assay_description  TEXT,
  assay_organism     TEXT,
  assay_format       TEXT,            -- ChEMBL BAO format, e.g. organism-based, microsome, cell-based
  system             TEXT NOT NULL CHECK (system IN ('in vivo', 'in vitro', 'unspecified')),
  -- Study details exactly as ChEMBL reports them (value first, then assay); blank when not given.
  dosed_drug         TEXT,            -- drug actually given, e.g. a salt, prodrug or the parent of a metabolite
  interacting_drug   TEXT,            -- co-administered drug
  health_status      TEXT,            -- e.g. healthy, cancer
  gender             TEXT,
  n_subjects         TEXT,
  age_range          TEXT,
  dietary_status     TEXT,            -- fed / fasting
  dose               TEXT,
  route              TEXT,
  regimen            TEXT,            -- single / multiple dose
  formulation        TEXT,
  sample_tissue      TEXT,            -- plasma, urine, bile
  document_chembl_id TEXT,
  PRIMARY KEY (compound_id, activity_id)
);

-- Binding records (assay type B) that have a pChEMBL value.
CREATE TABLE chembl_binding (
  compound_id        INTEGER NOT NULL REFERENCES compounds (compound_id),
  activity_id        INTEGER NOT NULL,
  molecule_chembl_id TEXT NOT NULL,
  assay_chembl_id    TEXT,
  assay_description  TEXT,
  standard_type      TEXT,
  standard_relation  TEXT,
  standard_value     REAL,
  standard_units     TEXT,
  pchembl_value      REAL,
  target_chembl_id   TEXT,
  target_pref_name   TEXT,
  target_organism    TEXT,
  document_chembl_id TEXT,
  PRIMARY KEY (compound_id, activity_id)
);

CREATE VIEW chembl_target_summary AS
  SELECT compound_id, target_chembl_id, target_pref_name, target_organism,
         COUNT(*) AS n_activities, AVG(pchembl_value) AS mean_pchembl, MAX(pchembl_value) AS max_pchembl
  FROM chembl_binding
  GROUP BY compound_id, target_chembl_id;

-- ice.py: plasma fraction unbound used by the ICE IVIVE tool (httk measured, OPERA predicted).
CREATE TABLE ice_fraction_unbound (
  ice_id         INTEGER PRIMARY KEY,
  compound_id    INTEGER NOT NULL REFERENCES compounds (compound_id),
  assay          TEXT NOT NULL,   -- e.g. 'httk, Human Plasma Fraction Unbound'
  endpoint       TEXT NOT NULL,   -- Fu, or Applicability_Domain for OPERA
  value          REAL,
  unit           TEXT,
  species        TEXT,
  reference      TEXT,
  substance_name TEXT,
  dtxsid         TEXT,
  casrn          TEXT,
  ice_inchikey   TEXT NOT NULL,   -- qsarReadyId returned by ICE
  match_method   TEXT NOT NULL    -- inchikey or inchikey_connectivity
);
