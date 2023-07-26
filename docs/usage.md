# Usage


## GUI

To call the script from the Graphical User Interface the folder containing the input files must be opened in terminal
and the correct python environment activated. The GUI is then simply called by:

```bash
python PyDP4_GUI.py
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
PyDP4.py Candidate CandidateNMR
```

where Candidate is the sdf file containing 3D coordinates of the candidate
structure (without the extension).

DP5 will automatically detect whether an NMR description (see below) or raw NMR data has been provided by the user.
If raw NMR data has been provided (currently in the Bruker file format), NMR-AI will automatically process the NMR spectra
provided. In this case the NMR data should be organised into a single folder (passed to DP5) containing one or two
folders labelled Proton and (or) Carbon.

Alternatively:

```bash
PyDP4.py -w gmns Candidate CandidateNMR
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
PyDP4.py -s chloroform Candidate CandidateNMR
```

If solvent is not given, no solvent is used.

### With explicit diastereomer/other candidate structures

```bash
PyDP4.py Candidate1 Candidate2 Candidate3 ... CandidateNMR
```

The script does not attempt to generate diastereomers, simply carries out the
DP4 on the specified candidate structures.

Structures can also be added from InChI, Smiles and Smarts strings using the designated switches. For example:

```bash
PyDP4.py --InChI Candidates_inchis.inchi ... CandidateNMR
```
where Candidates_inchs.inchi is a text file with all the desired InChI strings on separate lines.


### Other command line switches

Script has several other switches, including switching the molecular mechanics and dft software etc.

```
  -m {t,m}, --mm {t,m}  Select molecular mechanics program, t for tinker or m
                        for macromodel, default is t

  -d {j,g,n,z,w}, --dft {j,g,n,z,w}
                        Select DFT program, j for Jaguar, g for Gaussian, n
                        for NWChem, z for Gaussian on ziggy, w for NWChem on
                        ziggy, default is z (jaguar is not yet implemented)

  --StepCount STEPCOUNT
                        Specify stereocentres for diastereomer generation

  -s SOLVENT, --solvent SOLVENT
                        Specify solvent to use for dft calculations

  -q QUEUE, --queue QUEUE
                        Specify queue for job submission on ziggy
			(default is s1)

  -t NTAUT, --ntaut NTAUT
                        Specify number of explicit tautomers per diastereomer
                        given in structure files, must be a multiple of
                        structure files

  -r, --rot5            Manually generate conformers for 5-memebered rings

  --ra RA               Specify ring atoms, for the ring to be rotated, useful
                        for molecules with several 5-membered rings

  --AssumeDFTDone       Assume RMSD pruning, DFT setup and DFT calculations
                        have been run already (saves time when repeating DP4
			analysis)

  -g, --GenOnly         Only generate diastereomers and tinker input files,
                        but don't run any calculations (useful for diastereomer
			generation for calculations ran on computers
			without OpenBabel)

  -c STEREOCENTRES, --StereoCentres STEREOCENTRES
                        Specify stereocentres for diastereomer generation

  -T, --GenTautomers    Automatically generate tautomers

  -o, --DFTOpt          Optimize geometries at DFT level before NMR prediction

  --pd                  Use python port of DP4

  -b BASICATOMS, --BasicAtoms BASICATOMS
                        Generate protonated states on the specified atoms and
                        consider as tautomers
```

More information on those can be obtained by running `PyDP4.py -h`
