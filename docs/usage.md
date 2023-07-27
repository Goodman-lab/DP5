# Usage

The tool can be run in two modes:

- With a GUI, using the `PyDP4_GUI.py` script
- From the command line, using the `PyDP4.py` script


### Ways to play!

!!! warning

    Following either the [installation guide](installation.md) or [quickstart instructions](quickstart.md),
    you should have already created and activated a virtual environment containing the dependencies for the tool.

    **The virtual environment must be activated to run the tool, and needs to be re-activated for each new session!**



Can be run via lots of different commands





## GUI

To call the script from the Graphical User Interface the folder containing the input files must be opened in terminal
and the correct python environment activated. The GUI is then simply called by:

```bash
pydp4_gui
```

This will open the main DP5 GUI window, the user can then select the required input files and the settings for MM and
DFT calculations. The expected structure file format is *.sdf

Once the DP5 has finished all calculations the GUI will open a number of new tabs. These tabs display interactive
carbon and proton NMR spectra as processed and assigned by NMR-AI as well as interactive plots of the NMR prediction
errors and probabilities and conformer energies and populations.

## Terminal

To call the script from terminal:

### With all diastereomer generation

```bash
pydp4 Candidate CandidateNMR
```

where Candidate is the sdf file containing 3D coordinates of the candidate
structure (without the extension).

DP5 will automatically detect whether an NMR description (see below) or raw NMR data has been provided by the user.
If raw NMR data has been provided (currently in the Bruker file format), NMR-AI will automatically process the NMR spectra
provided. In this case the NMR data should be organised into a single folder (passed to DP5) containing one or two
folders labelled Proton and (or) Carbon.

Alternatively:

```bash
pydp4 -w gmns Candidate CandidateNMR
```

The -w switch specifies the PyDP4 workflow, c specifies structure cleaning utilising RDkit, g specifies diastereomer
generation, m molecular mechanics conformational searching, o DFT geometry optimisation, n DFT NMR calculation,
e M062X single point energy calculations, s calculate DP4 statistics and w for DP5 statistics.
The default workflow is gnms, optimum workflow is gnomes. The letter a can be given in place of s, in the instance
DP5 will analyse and assign the provided NMR spectra but not calculate any probabilities.

In addition the -s switch can be used to specify the solvent for use in MM and DFT calculations. Supported solvents are
listed in the TMSdata file. Other Gaussian/Jaguar/NWChem supported solvents can be used, but only with manually interpreted
data.

```bash
pydp4 -s chloroform Candidate CandidateNMR
```

If solvent is not given, no solvent is used.

### With explicit diastereomer/other candidate structures

```bash
pydp4 Candidate1 Candidate2 Candidate3 ... CandidateNMR
```

The script does not attempt to generate diastereomers, simply carries out the
DP4 on the specified candidate structures.

Structures can also be added from InChI, Smiles and Smarts strings using the designated switches. For example:

```bash
pydp4 --InChI Candidates_inchis.inchi ... CandidateNMR
```
where Candidates_inchs.inchi is a text file with all the desired InChI strings on separate lines.


### Command line flags

A full description of the rest of the command line flags, which include
switching the molecular mechanics and dft software etc. can be found on
the tool's manual/help page with the following command:

```bash
pydp4 --help
```


