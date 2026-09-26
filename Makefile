.DEFAULT_GOAL := help
.PHONY: help setup benchmark demo api web check clean-results

help: ## Show available targets
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[1m%-12s\033[0m %s\n", $$1, $$2}'

setup: ## Install everything and fetch the MHCflurry models and Ensembl proteome
	uv sync --extra presentation --extra plots --extra benchmark --extra api
	uv run mhcflurry-downloads fetch models_class1_presentation
	uv run neoantigene fetch-proteome

benchmark: ## Build the Ott 2017 case from the publisher and score it label-blind
	uv run --extra benchmark python scripts/build_ott2017_benchmark.py
	uv run python scripts/run_ott2017_benchmark.py

demo: ## Bring up the API and the web UI together
	./scripts/dev.sh

api: ## Run the FastAPI wrapper alone
	uv run --extra api uvicorn neoantigene_api.main:app --reload --port 8000 \
		--app-dir apps/api

web: ## Run the Next.js UI alone (expects the API on :8000)
	cd apps/web && npm run dev

check: ## Lint, typecheck and test
	uv run ruff check src scripts
	uv run ruff format --check src scripts
	uv run mypy
	uv run pytest -q

clean-results: ## Remove generated run outputs
	rm -rf results/
