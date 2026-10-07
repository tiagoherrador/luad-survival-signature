DATA ?= data
.PHONY: install test run
install:
	pip install -e ".[dev]"
test:
	pytest -q
run:
	luad-run --data-dir $(DATA) --out results --n-splits 30
