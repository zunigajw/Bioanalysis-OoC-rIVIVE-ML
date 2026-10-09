# Database Curation: `just` builds the Podman container from the Containerfile and runs GenDB inside it.
# Output: output/ooc_drug_ml.db and one spreadsheet per database (rdkit, chembl, ice).
default: database_curation

database_curation:
    podman run --rm -v "$PWD":/work -w /work $(podman build -q -f Containerfile .) make -f GenDB

clean:
    podman run --rm -v "$PWD":/work -w /work $(podman build -q -f Containerfile .) make -f GenDB clean
