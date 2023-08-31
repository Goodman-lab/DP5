# DP5

## What is DP5?

DP5 as described in our [latest paper](https://doi.org/10.1039/D1SC04406K) calculates the standalone
probability of structures being correct. If a user cannot guarantee the correct structure is in the list of proposals,
DP5 calculates the probability each is correct, quantifying this uncertainty.

[DP5](https://doi.org/10.1039/D1SC04406K) integrates NMR-AI, software for automatic processing, assignment
and visualisation of raw NMR data. This functionality affords fully automated DP4 analysis of databases of molecules
with no user input. DP5 also utilises a fully revised version of PyDP4 updated for Python 3 and with major workflow
improvements.

Users can now also choose to run DP5 within a Graphical user interface for increased ease of use when single DP4
calculations. This graphical user interface can also be used to explore the spectra processed and assigned by NMR-AI.
DP5 can also be utilised to perform assignments of NMR spectra for a single isomer of a molecule. These assignments
can be viewed in the graphical output from DP5 or investigated interactively using The included GUI application.

More details about the software can be found in the publication
[DP5 Straight from Spectrometer to Structure](https://doi.org/10.1039/D1SC04406K) and in the accompanying
supporting information.

## Documentation

[Documentation for this project is available as a website.](https://www.edmundgoodman.co.uk/DP5/)

TODO: Update link when merged into main repository!


## Quickstart

The following code can be copied and pasted into any bash terminal,
and should quickstart you running the tool. See the documentation for
more details, and specific system requirements.

```bash

# Download the source code from GitHub, and navigate into its directory
git clone https://github.com/Goodman-lab/DP5
cd DP5

# OPTIONAL; If this change is still in the restructure branch not master
git fetch origin
git checkout remotes/origin/restructure
# END OPTIONAL

# Run the make installation target
make install

# Activate the virtual environment containing the installed dependencies.
# Poetry is a python tool like pip which is installed by the make target if
# if isn't installed already
poetry shell

# Run the tool to see its manual. More instructions can then be found on the
# usage page on this website
pydp4 --help
```

## About the authors

DP5 was created by Alexander Howarth, Kristaps Ermanis, and Jonathan M. Goodman as part of the
[Goodman Lab](https://github.com/Goodman-lab/).
