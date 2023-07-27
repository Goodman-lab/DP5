# Implementation

## Code organisation

- `PyDP4.py` is the file that should be called to start the DP5 workflow. Interprets the
  arguments and takes care of the general workflow logic.

- `PyDP4_GUI.py` is the file that should be called to start the DP5 graphical user interface

- `InchiGen.py` gets called if diastereomer and/or tautomer and/or protomer generation is
  used. Called by `PyDP4.py`.

- `FiveConf.py` gets called if automatic 5-membered cycle corner-flipping is used. Called by
  `PyDP4.py`.

- `MacroModel.py` contains all of the MacroModel specific code for input generation, calculation
  execution and output interpretation. Called by `PyDP4.py`.

- `Tinker.py` contains all of the Tinker specific code for input generation, calculation
  execution and output interpretation. Called by `PyDP4.py`.

- `ConfPrune.pyx` is the Cython file for conformer alignment and RMSD pruning. Called by `Gaussian.py`
  and `NWChem.py`

- `Gaussian.py` contains all of the Gaussian specific code for input generation and calculation
  execution. Called by `PyDP4.py`.

- `GaussianZiggy.py` contains code specific to running Gaussian on Cambridge CMI local cluster. Called by `PyDP4.py`.

- `GaussianDarwin.py` contains code specific to running Gaussian on Cambridge CSD3 HPC cluster. Called by `PyDP4.py`.

- `NWChem.py` contains all of the NWChem specific code for input generation and calculation
  execution. Called by PyDP4.py.

- `NMRDP4GTF.py` takes care of all the NMR description interpretation, equivalent atom
  averaging, Boltzmann averaging, tautomer population optimisation (if used)
  and DP4 input preparation and running either `DP4.jar` or `DP4.py`. Called by
 `PyDP4.py`

- `nmrPredictNWChem.py` extracts NMR shifts from NWChem output files

- `DP4.py` is an equivalent and compact port of the original DP4 java implementation to Python
  of the same DP4 process. Has been updated to include multigaussian models.

- `DP5.py` is a python implementation of the DP5 probability calucaltion

- `Proton_processing.py` and `Carbon_processing.py` are NMR-AI python scripts from processing of raw proton/carbon NMR data

- `Proton_assignment.py` and `Carbon_assignment.py` are NMR-AI python scripts for assignment of raw proton/carbon NMR data

- `Proton_plotting.py` and `Carbon_plotting` are NMR-AI python script for plotting of raw proton/carbon NMR data processed and assigned by NMR-AI

- `StructureInput.py` is a script for cleaning structures using RDkit, generating 3d geometries and reading InChI, Smiles and SMARTS input formats
