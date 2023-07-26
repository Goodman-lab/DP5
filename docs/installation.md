# Installation



!!! note TODO

    Update this!

## Old instructions

All the python  files and one utility to convert to and from TINKER
nonstandard xyz file are in the attached archive. They are set up to work
from a centralised location. It is not a requirement, but it is probably
best to add the location to the PATH variable.

The script currently is set up to use MacroModel for molecular mechanics and
Gaussian for DFT and it runs Gaussian on ziggy by default. NWChem and TINKER is also supported.
This setup has several requirements.

1. One should have MacroModel or TINKER and NWChem or Gaussian.
   The DP5 folder can contain `settings.cfg`, where the locations for the installed software packages
   can be given, so that DP5 can launch them and use in the workflow. Examples on how to do
   this are provided in the included `settings.example` - this should be renamed to `settings.cfg` and
   any relevant changes made.

2. In order to run DP5 a python 3.6+ environment the following packages are required:
   - numpy
   - scipy
   - Cython
   - matplotlib
   - nmrglue
   - statsmodels
   - lmfit
   - openbabel
   - rdkit
   - pathos
   - qml
   - PyQT5 (optional if not using GUI)

3. RDKit and OpenBabel are required for the automatic diastereomer generation and
   other manipulations of sdf files (renumbering, ring corner flipping), as well as molecule
   graphical plotting. For OpenBabel this requirement includes Python bindings.
   The following links provide instructions for building OpenBabel with Python bindings:
   - http://openbabel.org/docs/dev/UseTheLibrary/PythonInstall.html
   - http://openbabel.org/docs/dev/Installation/install.html#compile-bindings
   Alternatively, the OpenBabel package in conda repository has also lately become
   more reliable, and is much easier to install via `conda install openbabel -c conda-forge`.
   Both OpenBabel 2.4.1 and 3.x are supported.

4. Finally, to run calculations on a computational cluster, a passwordless
   ssh connection should be set up in both directions:
   desktop -> cluster and cluster -> desktop. In most cases the modification
   of the relevant functions in `Gaussian.py` or `NWChem.py` will be required
   to fit your situation.
 
5. All development and testing was done on Linux. However, both the scripts
   and all the backend software should work equally well on MacOS with minor adjustments.
   Windows will require more modification, and is currently not supported.
