# mcycle-pipeline make targets.
#
# `make smoke` is the end-to-end test that exists today: it runs the pipeline on the
# 18-genome smoke panel (../test_panel) and checks the diagnostic calls against
# validation/smoke_expectations.tsv. It shows the pipeline runs and still gets the
# textbook organisms right. It is NOT an accuracy estimate — the accuracy campaign
# (reference panel + ground truth + hold-out + regression floors, as in the nitrogen
# and sulfur sister pipelines) is still to be done; see ROADMAP.md.

PANEL := ../test_panel

.PHONY: env panel ref-panel validate-panel dbs smoke smoke-check clean clean_all

env:
	@mamba env create -f envs/mcycle.yaml 2>/dev/null || conda env create -f envs/mcycle.yaml

# Download / copy the smoke-panel genomes listed in ../test_panel/panel.tsv.
panel:
	$(PANEL)/fetch_panel.sh

# Download the reference panel (validation/panel.tsv) into ../ref_panel and check that
# every file is the organism and the assembly the roster says it is.
ref-panel:
	validation/fetch_ref_panel.sh

validate-panel:
	python validation/validate_panel.py --online --strict

# (Re)build the HMM and BLAST databases from config/targets.yaml.
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
