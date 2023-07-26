### Variables ###
# https://clarkgrubb.com/makefile-style-guide
MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.DEFAULT_GOAL := install


### Installation ###
.PHONY: requirements.txt
requirements.txt:
	@echo "Exporting requirements.txt file"
	@poetry lock
	@poetry export \
		--without-hashes --format=requirements.txt \
		> requirements.txt

.PHONY: install
install:
	@echo "Initialising project:"
	@python3 -m venv .venv
	@.venv/bin/pip install -r requirements.txt
