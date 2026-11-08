# Pinned BLAST seeds

`unstable_refs.fasta` (19 records) and `blast_gated_refs.fasta` (9 records): the 19
curated UniProt accessions of `config/targets.yaml` (mcrA, mcrA_anme, pmoA, mmoX; pmoA and
mcrA_anme sit in both files), as fetched on 2026-11-07. Headers are
`>{target_id}||{accession}` followed by the UniProt header. mcycle never fetches seeds by
gene name, so every record here is one the curator listed.

**Which fetch this is.** The databases of the validation campaign (2026-10-05) were built
from a fetch of 2026-10-04; those files were overwritten by the rebuild of 2026-11-07
(review phase R6, `WORKPLAN_2026-11_review.md`). The two fetches differ in the UniProt
description text of the headers only: the 36 diamond tables of the 18-genome smoke panel
and all 54 call tables are byte-identical between the two databases, and the 49-genome
reference panel passes its 14 regression checks with unchanged score tables. The
sequences are the validated ones; the header text is UniProt's of 2026-11-07.

**Why they are in the repository.** UniProt entries are revised and deleted, so a
database rebuilt from UniProt is not guaranteed to be the validated one, and a rebuild
(triggered by any edit of `config/targets.yaml`) needed the network.
`workflow/scripts/build_blast_db.py` copies these files into `resources/blast_db/` and
indexes them. With `--upstream`, or when `config/targets.yaml` lists a curated accession
that is not in these files, it fetches all seeds from UniProt instead — after adding a
seed, rebuild that way, check the result (`make smoke`, `make regression`), and replace
the two files here.

Sequences are from UniProtKB (CC BY 4.0).
