### Variables ###
# https://clarkgrubb.com/makefile-style-guide
MAKEFLAGS += --warn-undefined-variables
SHELL := bash
.DEFAULT_GOAL := install


### Installation ###
wheels:
	@echo "Downloading wheels"
	@pip wheel --use-pep517 "qml (==0.4.0.27)" &&\
		mkdir -p lib/ &&\
		mv ./qml-0.4.0.27-*.whl wheels/

.PHONY: install
install: wheels
	@echo "Initialising project:"
	@pip install poetry
	@poetry install
