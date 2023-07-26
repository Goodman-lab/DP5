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

