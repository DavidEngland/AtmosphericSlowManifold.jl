# ==============================================================================
# Makefile for AtmosphericSlowManifold.jl
# ==============================================================================

JULIA       ?= julia
PROJECT     ?= --project=.
SMOKE_ENV   ?= ASM_RUN_SMOKE=1

.PHONY: help instantiate update test test-smoke test-calibration test-export e2e campaign-export campaign-summary campaign-validate campaign-clean pde-benchmark report clean

default: help

## help: Display this help menu with available targets
help:
	@echo "AtmosphericSlowManifold.jl - Development Workflow Targets:"
	@echo ""
	@sed -n 's/^##//p' $(MAKEFILE_LIST) | column -t -s ':' | sed -e 's/^/ /'
	@echo ""

## instantiate: Resolve and instantiate Julia project dependencies
instantiate:
	$(JULIA) $(PROJECT) -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'

## update: Update project dependencies
update:
	$(JULIA) $(PROJECT) -e 'using Pkg; Pkg.update()'

## test: Run the full standard test suite
test:
	$(JULIA) $(PROJECT) -e 'include("test/runtests.jl")'

## test-smoke: Run full test suite including solver-level smoke tests (ASM_RUN_SMOKE=1)
test-smoke:
	$(SMOKE_ENV) $(JULIA) $(PROJECT) -e 'include("test/runtests.jl")'

## test-calibration: Run calibration interface & backend tests
test-calibration:
	$(JULIA) $(PROJECT) -e 'include("test/test_calibration_interface.jl")'

## test-export: Run multi-format export utilities test suite
test-export:
	$(JULIA) $(PROJECT) -e 'include("test/test_export_utilities.jl")'

## e2e: Run end-to-end pipeline test (Ingestion -> Discovery -> PDE -> Calibration -> Export)
e2e:
	$(JULIA) $(PROJECT) -e 'include("test/test_end_to_end_pipeline.jl")'

## campaign-export: Run end-to-end multi-campaign data extraction & export pipeline
campaign-export:
	$(JULIA) $(PROJECT) scripts/run_campaign_exports.jl

## campaign-summary: Display generated campaign summary report table
campaign-summary:
	@cat reports/generated/campaign_exports/tables/campaign_summary.md

## campaign-validate: Validate generated campaign artifacts, schemas, and report sections
campaign-validate:
	$(JULIA) $(PROJECT) scripts/validate_campaign_exports.jl

## campaign-clean: Remove generated campaign export artifacts
campaign-clean:
	rm -rf reports/generated/campaign_exports/

## pde-benchmark: Run physical-closure spectral PDE benchmark and export summary artifacts
pde-benchmark:
	$(JULIA) $(PROJECT) scripts/run_pde_closure_benchmark.jl

## report: Build benchmark PDF report from template and generated artifacts
report:
	@mkdir -p reports/generated
	cp templates/report.tex.mustache reports/generated/report.tex
	cd reports/generated && pdflatex -interaction=nonstopmode report.tex > /dev/null 2>&1
	cd reports/generated && pdflatex -interaction=nonstopmode report.tex > /dev/null 2>&1
	@echo "PDF build complete: reports/generated/report.pdf"

## clean: Remove generated output artifacts, figures, and temporary files
clean:
	rm -rf reports/
	rm -f *.nc *.csv *.json
	find . -type f -name "*.cov" -delete
	find . -type f -name "*.mem" -delete