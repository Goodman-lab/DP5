### Variables ###
# https://clarkgrubb.com/makefile-style-guide
MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.DEFAULT_GOAL := install


### Installation ###
poetry.lock:
	@poetry lock

requirements.txt: poetry.lock
	@poetry export \
		--without-hashes --format=requirements.txt \
		> requirements.txt

.PHONY: install
install: requirements.txt
	@echo "Initialising project:"
	@python3 -m venv .venv
	@.venv/bin/pip install -r requirements.txt
