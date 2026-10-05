# mcycle-pipeline make targets.
#
# `make smoke` runs the pipeline on the 18-genome smoke panel (test_panel/) and checks
# the diagnostic calls against validation/smoke_expectations.tsv. It shows the pipeline
# runs and still gets the textbook organisms right. It is NOT an accuracy estimate —
# that is `make regression` (49-genome reference panel, floors of
# validation/test_regression.py).
# Run the targets inside the conda env (`conda activate cycle-pipeline`).

PANEL := test_panel

.PHONY: env panel ref-panel validate-panel ground-truth regression regression-score dbs smoke smoke-check clean clean_all

env:
	@mamba env create -f envs/mcycle.yaml 2>/dev/null || conda env create -f envs/mcycle.yaml

# Download / copy the smoke-panel genomes listed in test_panel/panel.tsv.
panel:
	$(PANEL)/fetch_panel.sh

# Download the reference panel (validation/panel.tsv) into ../ref_panel and check that
# every file is the organism and the assembly the roster says it is.
ref-panel:
	validation/fetch_ref_panel.sh

validate-panel:
	python validation/validate_panel.py --online --strict

# Rebuild both ground truths: KEGG-derived, then curated-function (KEGG + the cells of
# validation/curated_cells.tsv + mcrA_anme from validation/phenotype_gt.tsv).
ground-truth:
	python validation/build_ground_truth.py
	python validation/build_curated_function_gt.py

# Accuracy regression gate: run the 49-genome reference panel (own config and results
# directory), re-check seed leakage, score, and fail below the floors of
# validation/test_regression.py. `regression-score` re-scores existing results_ref/.
regression: ref-panel
	MCYCLE_CONFIG=config/config_ref.yaml python run.py --input ../ref_panel --prodigal-mode single --cores 8 --skip-db-setup
	$(MAKE) regression-score

regression-score:
	python validation/detect_seed_leakage.py
	python validation/score_mcycle.py
	python validation/test_regression.py

# (Re)build the HMM and BLAST databases from config/targets.yaml (KOfam profiles from
# resources/kofam_pinned/, seeds from UniProt).
dbs:
	python workflow/scripts/build_hmm_db.py --force
	python workflow/scripts/build_blast_db.py --force

smoke: panel
	python run.py --input $(PANEL) --prodigal-mode single --cores 8
	python validation/check_smoke.py

# Re-check the EXISTING results/ without running the pipeline.
smoke-check:
	python validation/check_smoke.py

clean:
	rm -rf results/*

clean_all: clean
	rm -rf resources/hmm/* resources/blast_db/*
