### Variables ###
# https://clarkgrubb.com/makefile-style-guide
MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.DEFAULT_GOAL := install


### Installation ###
.PHONY: install
install: poetry
	@echo "Initialising project:"
	@poetry install --without=qml
	@.venv/bin/pip install qml  # Install QML separately to avoid build dependency issues

# Later can activate pre-commit on install
#	@poetry run pre-commit install

.PHONY: requirements.txt
requirements.txt: poetry
	@echo "Exporting requirements.txt file"
	@poetry lock
	@poetry export \
		--without-hashes --format=requirements.txt \
		> requirements.txt

.PHONY: poetry
poetry:
	@echo "Ensuring poetry is installed"
	@poetry --version || pip install poetry


### Static analysis ###
.PHONY: check
check: .venv/
	@echo "Running pre-commit hooks:"
	@poetry run pre-commit run --all-files


### Documentation ###
.PHONY: docs
docs: .venv/
	@echo "Generating documentation:"
	@poetry run mkdocs build --strict


### Cleanup ###
.PHONY: clean_docs
clean_docs:
	@echo "Deleting generated documentation site:"
	rm -rf site/

.PHONY: clean_pycache
clean_pycache:
	@echo "Deleting python cache files:"
	find . -not -path "./.venv/*" | \
		grep -E "(/__pycache__$$|\.pyc$$|\.pyo$$)" | \
		xargs rm -rf

.PHONY: clean
clean: clean_docs clean_pycache

.PHONY: clobber
clobber: clean
	@echo "Deleting virtual environment"
	rm -rf .venv/
