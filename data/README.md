# data/ — raw inputs (read-only)

Raw data in this directory is **immutable and read-only** by convention; a hook blocks writes here.

Place the dataset and its metadata file here. Nothing in this directory is ever modified, cleaned, or overwritten in place.

Everything in `results/` and `figures/` is regenerated from the contents of this directory plus a script — derived outputs never live here.
