#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
PyDP4 integrated workflow for the running of MM, DFT GIAO calculations and
DP4 analysis
v1.0

Copyright (c) 2015-2019 Kristaps Ermanis, Alexander Howarth, Jonathan M. Goodman

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Created on Wed Nov 19 15:26:32 2014
Updated on Feb 07 2019

@author: ke291

The main file, that should be called to start the PyDP4 workflow.
Interprets the arguments and takes care of the general workflow logic.
"""

from dp5 import NMR
from dp5 import DP5
from dp5 import DP4

import sys
import os
import datetime
import argparse
import importlib
import getpass
from pathlib import Path

DFTpackages = [['n', 'w', 'g', 'z', 'd'],['NWChem', 'NWChemZiggy', 'Gaussian', 'GaussianZiggy', 'GaussianDarwin']]
MMpackages = [['m','t'],['Macromodel','Tinker']

if os.name == 'nt':
    import pyximport

    pyximport.install()
    # import dp5.ConfPrune as ConfPrune
    from dp5 import ConfPrune
else:
    import pyximport

    pyximport.install()
    # import dp5.ConfPrune as ConfPrune
    from dp5 import ConfPrune

# Assigning the config default values
# Settings are defined roughly in the order they are used in the script
class Settings:
    # --- Main options ---
    MM = 'm'  # m for MacroModel, t for Tinker
    DFT = 'z'  # n, g, z or for NWChem or Gaussian
    Workflow = 'gmns'  # defines which steps to include in the workflow
    # c for RDkit cleaning of input structures and 3d coordinate generation
    # g for generate diastereomers
    # m for molecular mechanics conformational search
    # o for DFT optimization
    # e for DFT single-point energies
    # n for DFT NMR calculation
    # s for computational and experimental NMR data extraction and stats analysis
    # w for DP5 probability calculation
    Solvent = ''  # solvent for DFT optimization and NMR calculation
    ScriptDir = ''  # Script directory, automatically set on launch
    InputFiles = []  # Structure input files - can be MacroModel *-out.mae or *sdf files
    InputFilesPaths = []  # Path object for Structure input files - can be MacroModel *-out.mae or *sdf files
    NMRsource = ''  # File or folder containing NMR description or data
    Title = 'DP4molecule'  # Title of the calculation, set to NMR file name by default on launch
    AssumeDone = False  # Assume all computations done, only read DFT output data and analyze (use for reruns)
    AssumeConverged = False  # Assume all optimizations have converged, do NMR and/or energy calcs on existing DFT geometries
    UseExistingInputs = False  # Don't regenerate DFT inputs, use existing ones. Good for restarting a failed calc
    Smiles = None  # Smiles input file - text file with smiles strings on separate lines
    InChIs = None  # InChI input file - text file with inchi strings on separate lines
    Smarts = None  # Smarts input file - text file with Smarts strings on separate lines

    # --- Diastereomer generation ---
    SelectedStereocentres = []  # which stereocentres to vary for diastereomer generation

    # --- Molecular mechanics ---
    ForceField = 'mmff'  # ff tfOPto use for conformational search
    MMstepcount = 10000  # Max number of MM steps to do, if less than MMfactor*rotable_bonds
    MMfactor = 2500  # MMfactor*rotable_bonds gives number of steps to do if less than MMstepcount
    Rot5Cycle = False  # Special dealing with 5-membered saturated rings, see FiveConf.py
    RingAtoms = []  # Define the 5-membered ring, useful if several are present in molecule
    SCHRODINGER = ''  # Define the root folder for Schrodinger software
    TinkerPath = '/tinker'  # Define the root folder for Tinker software,
    # must contain bin/scan and params/mmff.prm for the process to work

    # --- Conformer pruning ---
    HardConfLimit = 1000  # Immediately stop if conformers to run exceed this number
    ConfPrune = True  # Should we prune conformations?
    PerStructConfLimit = 100  # Max numbers of conformers allowed per structure for DFT stages
    InitialRMSDcutoff = 0.75  # Initial RMSD threshold for pruning
    MaxCutoffEnergy = 10.0  # Max conformer MM energy in kJ/mol to allow

    # --- DFT ---
    NWChemPath = "nwchem"  # Path to nwchem executable. If it's in the path, can be just 'nwchem'
    GausPath = ""  # Path to Gaussian executable. If it's in the path, can be just 'g09' or 'g16'
    # If left empty, it will attempt to use g09 in GAUS_EXEDIR environment variable
    MaxDFTOptCycles = 50  # Max number of DFT geometry optimization cycles to request.
    CalcFC = False  # Calculate QM force constants before optimization
    OptStepSize = 30  # Max step Gaussian should take in geometry optimization
    charge = None  # Manually specify charge for DFT calcs
    nBasisSet = "6-311g(d)"  # Basis set for NMR calcs
    nFunctional = "mPW1PW91"  # Functional for NMR calcs
    oBasisSet = "6-31g(d,p)"  # Basis set for geometry optimizations
    oFunctional = "b3lyp"  # Functional for geometry optimizations
    eBasisSet = "def2tzvp"  # Basis set for energy calculations
    eFunctional = "m062x"  # Functional for energy calculations

    # --- Computational clusters ---
    """ These should probably be moved to relevant *.py files as Cambridge specific """
    user = ''  # Linux user on computational clusters, not used for local calcs
    TimeLimit = 24  # Queue time limit on comp clusters
    queue = 'SWAN'  # Which queue to use on Ziggy
    project = 'GOODMAN-SL3-CPU'  # Which project to use on Darwin
    DarwinScrDir = '/home/u/rds/hpc-work/'  # Which scratch directory to use on Darwin
    StartTime = ''  # Automatically set on launch, used for folder names
    nProc = 1  # Cores used per job, must be less than node size on cluster
    DarwinNodeSize = 32  # Node size on current CSD3
    MaxConcurrentJobsZiggy = 75  # Max concurrent jobs to submit on ziggy
    MaxConcurrentJobsDarwin = 320  # Max concurrent jobs to submit on CSD3

    # --- NMR analysis ---
    TMS_SC_C13 = 191.69255  # Default TMS reference C shielding constant (from B3LYP/6-31g**)
    TMS_SC_H1 = 31.7518583  # Default TMS reference H shielding constant (from B3LYP/6-31g**)

    # --- Stats ---
    StatsModel = 'g'  # What statistical model type to use
    StatsParamFile = 'none'  # Where to find statistical model parameters

    # --- Output folder ---
    OutputFolder = ''  # folder to print dp4 output to - default is cwd
    GUIRunning = False  # Boolean has PyDP4 been called from commandline or from GUI


settings = Settings()


# Data structure keeping all of isomer data in one place.
class Isomer:
    def __init__(self, InputFile, Charge=-100):
        self.InputFile = InputFile  # Initial structure input file
        self.BaseName = InputFile  # Basename for other files
        self.Atoms = []  # Element labels
        self.Conformers = []  # from conformational search, list of atom coordinate lists
        self.MMCharge = 0  # charge from conformational search
        self.ExtCharge = Charge  # externally provided charge
        self.RMSDCutoff = 0  # RMSD cutoff eventually used to get the conformer number below the limit
        self.DFTConformers = []  # from DFT optimizations, list of atom coordinate lists
        self.ConfIndices = []  # List of conformer indices from the original conformational search for reference
        self.MMEnergies = []  # Corresponding MM energies in kj/mol
        self.DFTEnergies = []  # Corresponding DFT energies in hartrees
        self.Energies = []  # Final energies used in conformer population prediction in kj/mol
        self.Populations = []  # Conformer populations
        self.OptInputFiles = []  # list of DFT NMR input file names
        self.OptOutputFiles = []  # list of DFT NMR output file names
        self.EInputFiles = []  # list of DFT NMR input file names
        self.EOutputFiles = []  # list of DFT NMR output file names
        self.NMRInputFiles = []  # list of DFT NMR input file names
        self.NMROutputFiles = []  # list of DFT NMR output file names
        self.ShieldingLabels = []  # A list of atom labels corresponding to the shielding values
        self.ConformerShieldings = []  # list of calculated NMR shielding constant lists for every conformer
        self.ConformerCShifts = [] # list of calculated C NMR shifts lists for every conformer
        self.ConformerHShifts = [] # list of calculated H NMR shifts lists for every conformer
        self.BoltzmannShieldings = []  # Boltzmann weighted NMR shielding constant list for the isomer
        self.Cshifts = []  # Calculated C NMR shifts
        self.Hshifts = []  # Calculated H NMR
        self.Clabels = []
        self.Hlabels = []
        self.Cexp = []  # Experimental C NMR shifts
        self.Hexp = []  # Experimental H NMR shifts

def main(settings):


    print("Current working directory: " + os.getcwd())
    print("Initial input files: " + str(settings.InputFiles))
    print("NMR file: " + str(settings.NMRsource))
    print("Workflow: " + str(settings.Workflow))

    # Read in any text inputs and add these to the input file list

    from dp5 import StructureInput

    if settings.Smiles:
        settings.InputFiles.extend(StructureInput.GenerateSDFFromTxt(settings.Smiles, 'Smiles'))

    if settings.Smarts:
        settings.InputFiles.extend(StructureInput.GenerateSDFFromTxt(settings.Smarts, 'Smarts'))

    if settings.InChIs:
        settings.InputFiles.extend(StructureInput.GenerateSDFFromTxt(settings.InChIs, 'InChI'))

    # Clean up input files if c in workflow - this generates a new set of 3d coordinates as a starting point

    if 'c' in settings.Workflow and len(settings.InputFiles) > 0:
        from dp5 import StructureInput

        # if requested generate 3d coordinates to define any stereochemistry

        settings.InputFiles = StructureInput.CleanUp(settings.InputFiles)

    # if no structure inputs have been found at this point quit

    if len(settings.InputFiles) == 0:
        print("\nNo input files were found please use -h for help with input options quitting...")

        quit()

    # if g in workflow check number of stereocentres for each input and generate and diastereomers

    if ('g' in settings.Workflow):

        from dp5 import InchiGen
        print("\nGenerating diastereomers...")

        FinalInputFiles = []

        nStereo = [StructureInput.NumberofStereoCentres(InpFile) for InpFile in settings.InputFiles]

        if len(settings.InputFiles) == 1:

            FinalInputFiles.extend(
                InchiGen.GenDiastereomers(settings.InputFiles[0], nStereo[0], settings.SelectedStereocentres))

        else:

            for InpFile, nStereoCentres in zip(settings.InputFiles, nStereo):
                FinalInputFiles.extend(InchiGen.GenDiastereomers(InpFile, nStereoCentres, []))

        settings.InputFiles = list(FinalInputFiles)

    settings.InputFilesPaths = [Path.cwd() / i for i in settings.InputFiles]

    print("Generated input files: " + str(settings.InputFiles) + '\n')

    # Create isomer data structures
    Isomers = [Isomer(f.split('.sdf')[0]) for f in settings.InputFiles]

    print("Assuming all computations are done? ... ", settings.AssumeDone)
    print("Using preexisting DFT data? ... ", settings.UseExistingInputs)

    # Run conformational search, if requested
    if ('m' in settings.Workflow) and not (settings.AssumeDone or settings.UseExistingInputs):
        MM = ImportMM(settings.MM)
        print("Performing conformational search")
        print('\nSetting up files...')
        MMInputs = MM.SetupMM(settings)
        print('MM inputs: '+str(MMInputs))
        print('Running MM:')
        MMOutputs = MM.RunMM(MMInputs)
        print('\nReading conformers...')
        Isomers = MM.ReadConformers(MMOutputs, Isomers, settings)
        for iso in Isomers:
            print(iso.InputFile + ": " + str(len(iso.Conformers)) + ' conformers read within energy window')
    else:
        print('No conformational search was requested. Skipping...')
        settings.ConfPrune = False

    # Prune conformations, if requested.
    # For each isomer, the conformers list is replaced with a smaller list of conformers
    if (settings.ConfPrune) and not (settings.AssumeDone or settings.UseExistingInputs):
        print('\nPruning conformers...')
        Isomers = ConfPrune.RMSDPrune(Isomers, settings)
        for iso in Isomers:
            print(iso.InputFile + ": " + str(len(iso.Conformers)) + ' conformers after pruning with ' +
                  str(iso.RMSDCutoff) + 'A RMSD cutoff')

    if ('n' in settings.Workflow) or ('o' in settings.Workflow) \
            or ('e' in settings.Workflow) or settings.AssumeDone:
        DFT = ImportDFT(settings.DFT)
    else:
        print('\nNo DFT calculations were requested. Skipping...')

    if not (settings.AssumeDone):

        # Run DFT optimizations, if requested
        if ('o' in settings.Workflow):

            now = datetime.datetime.now()
            settings.StartTime = now.strftime('%d%b%H%M')

            print('\nSetting up geometry optimization calculations...')
            Isomers = DFT.SetupOptCalcs(Isomers, settings)
            print('\nRunning geometry optimization calculations...')
            Isomers = DFT.RunOptCalcs(Isomers, settings)

            print('\nReading DFT optimized geometries...')

            Isomers = DFT.ReadGeometries(Isomers, settings)

            # Add convergence check here before continuing with calcs!
            if (DFT.Converged(Isomers) == False) and (settings.AssumeConverged == False):
                print('Some of the conformers did not converge, quitting...')
                quit()

        # Run DFT single-point energy calculations, if requested
        if ('e' in settings.Workflow):

            now = datetime.datetime.now()
            settings.StartTime = now.strftime('%d%b%H%M')

            print('\nSetting up energy calculations...')
            Isomers = DFT.SetupECalcs(Isomers, settings)
            print('\nRunning energy calculations...')
            Isomers = DFT.RunECalcs(Isomers, settings)
            print('\nReading data from the output files...')
            Isomers = DFT.ReadEnergies(Isomers, settings)
            print("Energies: ")
            for iso in Isomers:
                print(iso.InputFile + ": " + str(iso.DFTEnergies))

        # Run DFT NMR calculations, if requested
        if ('n' in settings.Workflow):

            now = datetime.datetime.now()
            settings.StartTime = now.strftime('%d%b%H%M')

            print('\nSetting up NMR calculations...')
            Isomers = DFT.SetupNMRCalcs(Isomers, settings)
            print('\nRunning NMR calculations...')
            Isomers = DFT.RunNMRCalcs(Isomers, settings)
            print('\nReading data from the output files...')
            Isomers = DFT.ReadShieldings(Isomers)
            print("Shieldings: ")
            for iso in Isomers:
                print(iso.InputFile + ": ")
                for conf in iso.ConformerShieldings:
                    print(str(conf))

            Isomers = DFT.ReadEnergies(Isomers, settings)
            print("Energies: ")
            for iso in Isomers:
                print(iso.InputFile + ": " + str(iso.DFTEnergies))

    else:
        # Read DFT optimized geometries, if requested
        if ('o' in settings.Workflow):
            Isomers = DFT.GetPrerunOptCalcs(Isomers)
        if ('e' in settings.Workflow):
            Isomers = DFT.GetPrerunECalcs(Isomers)
        if ('n' in settings.Workflow):
            Isomers = DFT.GetPrerunNMRCalcs(Isomers)

        Isomers = DFT.ReadGeometries(Isomers, settings)

        # Read DFT NMR data, if requested
        if ('n' in settings.Workflow):
            Isomers = DFT.ReadShieldings(Isomers)
            Isomers = DFT.ReadEnergies(Isomers, settings)

    if not (NMR.NMRDataValid(Isomers)) or ('n' not in settings.Workflow):
        print('\nNo NMR data calculated, quitting...')
        quit()

    if ('s' in settings.Workflow) or ('a' in settings.Workflow) or ('w' in settings.Workflow):

        print('\nSetting TMS computational NMR shielding constant references')
        settings.TMS_SC_C13, settings.TMS_SC_H1 = NMR.GetTMSConstants(settings)

        print('\nConverting DFT data to NMR shifts...')
        Isomers = NMR.CalcBoltzmannWeightedShieldings(Isomers)
        Isomers = NMR.CalcNMRShifts(Isomers, settings)

        print('\nReading experimental NMR data...')
        NMRData = NMR.NMRData(settings)

        """
        print("Conformation data:")
        NMR.PrintConformationData(Isomers)
        """

        if NMRData.Type == 'desc':

            print('Experimental NMR description found and read.')

            # performs a pairwise assignment

            Isomers = NMR.PairwiseAssignment(Isomers, NMRData)

            print('Cshifts: ' + str(NMRData.Cshifts))
            print('Hshifts: ' + str(NMRData.Hshifts))

            print('Equivalents: ' + str(NMRData.Equivalents))
            print('Omits: ' + str(NMRData.Omits))

        elif NMRData.Type == "fid":

            for f in settings.NMRsource:

                if f.name == "Proton" or f.name == "proton":

                    from dp5.Proton_assignment import AssignProton
                    from dp5.Proton_plotting import PlotProton

                    print('\nAssigning proton spectrum...')
                    Isomers = AssignProton(NMRData, Isomers, settings)

                    if settings.GUIRunning == False:
                        print('\nPlotting proton spectrum...')
                        PlotProton(NMRData, Isomers, settings)

                elif f.name == "Carbon" or f.name == "carbon":

                    from dp5.Carbon_assignment import AssignCarbon
                    from dp5.Carbon_plotting import PlotCarbon

                    print('\nAssigning carbon spectrum...')
                    Isomers = AssignCarbon(NMRData, Isomers, settings)

                    if settings.GUIRunning == False:
                        print('\nPlotting carbon spectrum...')
                        PlotCarbon(NMRData, Isomers, settings)

        elif NMRData.Type == "jcamp":

            for f in settings.NMRsource:

                if f.name == "Proton.dx" or f.name == "proton.dx":

                    from dp5.Proton_assignment import AssignProton
                    from dp5.Proton_plotting import PlotProton

                    print('\nAssigning proton spectrum...')
                    Isomers = AssignProton(NMRData, Isomers, settings)

                    if settings.GUIRunning == False:
                        print('\nPlotting proton spectrum...')
                        PlotProton(NMRData, Isomers, settings)

                elif f.name == "Carbon.dx" or f.name == "carbon.dx":

                    from dp5.Carbon_assignment import AssignCarbon
                    from dp5.Carbon_plotting import PlotCarbon

                    print('\nAssigning carbon spectrum...')
                    Isomers = AssignCarbon(NMRData, Isomers, settings)

                    if settings.GUIRunning == False:
                        print('\nPlotting carbon spectrum...')
                        PlotCarbon(NMRData, Isomers, settings)

            print('Raw FID NMR datafound and read.')

        # print('\nProcessing experimental NMR data...')

        # NMRdata = NMR.ProcessNMRData(Isomers, settings.NMRsource, settings)

    if 'w' in settings.Workflow:

        if "o" not in settings.Workflow:

            print( "DFT optimised geometries required for DP5 calculation, please rerun with o option in workflow...")

            quit()

        print('\nCalculating DP5 probabilities...')

        # make folder for WF data to go into

        DP5data = DP5.DP5data(Path(settings.ScriptDir), len(Isomers[0].Atoms))

        if not os.path.exists('dp5'):

            os.mkdir(Path(settings.OutputFolder) / 'dp5')

            DP5data = DP5.ProcessIsomers(DP5data, Isomers, settings)
            DP5data = DP5.InternalScaling(DP5data)
            DP5data = DP5.kde_probs(Isomers, DP5data, 0.025)
            DP5data = DP5.BoltzmannWeight_DP5(Isomers, DP5data)
            DP5data = DP5.Calculate_DP5(DP5data)
            DP5data = DP5.Rescale_DP5(DP5data, settings)
            DP5data = DP5.Pickle_res(DP5data, settings)

        else:

            DP5data = DP5.UnPickle_res(DP5data, settings)


        DP5data = DP5.MakeOutput(DP5data, Isomers, settings)

    else:

        DP5data = []

    if 's' in settings.Workflow:

        if len(Isomers) < 2:

            print("Multiple structures required for DP4 probability calculations...")

        else:

            print('\nCalculating DP4 probabilities...')

            DP4data = DP4.DP4data()
            DP4data = DP4.ProcessIsomers(DP4data, Isomers)
            DP4data = DP4.InternalScaling(DP4data)
            DP4data = DP4.CalcProbs(DP4data, settings)
            DP4data = DP4.CalcDP4(DP4data)

            DP4data = DP4.MakeOutput(DP4data, Isomers, settings)

    else:
        print('\nNo DP4 analysis requested.')

        DP4data = []

    print('\nPyDP4 process completed successfully.')

    print("workflow" , settings.Workflow)

    return NMRData, Isomers, settings, DP4data, DP5data

# Selects which MM package to import, returns imported module

def ImportMM(mm):
    if mm in MMpackages[0]:
        MMindex = MMpackages[0].index(mm)
        MM = importlib.import_module(f"dp5.{MMpackages[1][MMindex]}")
    else:
        print("Invalid MM package selected")
        quit()
    return MM

# Selects which DFT package to import, returns imported module
def ImportDFT(dft):
    if dft in DFTpackages[0]:
        DFTindex = DFTpackages[0].index(dft)
        DFT = importlib.import_module(f"dp5.{DFTpackages[1][DFTindex]}")
    else:
        print("Invalid DFT package selected")
        quit()

    return DFT


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def NMR_files(NMR_args):

    print("NMR_path")

    NMR_path = Path(NMR_args)
    NMR_Data = []

    # check if path is from cwd or elsewhere:
    if len(NMR_path.parts) == 1:
        # if so a folder in the cwd has been passed add the cwd to the path
        NMR_path = Path.cwd() / NMR_path
        print(NMR_path)

    # now check if it is a directory or a file, add proton and carbon data here
    if NMR_path.is_dir():

        p_switch = 0
        c_switch = 0

        for f in NMR_path.iterdir():

            if f.name == "Carbon" or f.name == "carbon" or f.name == "Carbon.dx" or f.name == "carbon.dx":
                NMR_Data.append(f)
                c_switch = 1

            elif f.name == "Proton" or f.name == "proton" or f.name == "Proton.dx" or f.name == "proton.dx":
                NMR_Data.append(f)
                p_switch = 1

            if p_switch == 1 and c_switch == 1:
                break

        # self.NMR_list.addItem(str(filename[0].split("/")[-1]))

        if p_switch == 0 and c_switch == 0:
            NMR_Data.append(f)

    # if its not a directory add the file

    else:

        NMR_Data.append(NMR_path)

    settings.NMRsource = NMR_Data

    return


# Read the config file and fill in the corresponding attributes in settings class
def ReadConfig(settings):
    cfgpath = os.path.join(getScriptPath(), 'settings.cfg')
    if not os.path.exists(cfgpath):
        print('settings.cfg is missing!')
        return settings

    configfile = open(cfgpath, 'r')
    config = configfile.readlines()
    configfile.close()

    # Read in the new settings values from config
    newsettings = []
    for line in config:
        if ('#' in line) or (len(line) < 3) or ('=' not in line):
            continue

        newsettings.append([x.strip() for x in line[:-1].split('=')])
        if len(newsettings[-1]) < 2:
            newsettings[-1].append('')

    # Set the attributes in the settings class
    print('Settings read from settings.cfg:')
    for setting in newsettings:
        if hasattr(settings, setting[0]):
            setattr(settings, setting[0], setting[1])
            print('  ' + setting[0] + ': ' + setting[1])

    return settings


def run():
    print("==========================")
    print("PyDP4 script,\nintegrating Tinker/MacroModel,")
    print("Gaussian/NWChem and DP4\nv1.0")
    print("\nCopyright (c) 2015-2019 Kristaps Ermanis, Alexander Howarth, Jonathan M. Goodman")
    print("Distributed under MIT license")
    print("==========================\n\n")

    # This is a big anti-pattern, but there is not a better way to provide an entrypoint function here without a
    # major refactor!!!
    global settings

    # Read config file and fill in settings in from that
    settings = ReadConfig(settings)

    # These are then overridden by any explicit parameters given through the command line
    parser = argparse.ArgumentParser(description='PyDP4 script to setup\
    and run Tinker, Gaussian (on ziggy) and DP4')
    parser.add_argument('-w', '--workflow', help="Defines which steps to include in the workflow, " +
                                                 "can contain g for generate diastereomers, m for molecular mechanics conformational search, " +
                                                 "o for DFT optimization, e for DFT single-point energies, n for DFT NMR calculation, " +
                                                 "a for computational and experimental NMR data extraction " +
                                                 "s for computational and experimental NMR data extraction and stats analysis, default is 'gmns'",
                        default=settings.Workflow)
    parser.add_argument('-m', '--mm', help="Select molecular mechanics program,\
    t for tinker or m for macromodel, default is m", choices=MMpackages[0],
                        default='m')
    parser.add_argument('-d', '--dft', help="Select DFT program, \
    g for Gaussian, n for NWChem, z for Gaussian on ziggy, d for Gaussian on \
    Darwin, default is g", choices=DFTpackages[0], default='g')
    parser.add_argument('--StepCount', help="Specify\
    stereocentres for diastereomer generation")

    parser.add_argument('StructureFiles', nargs='*', default=[], help=
    "One or more SDF file for the structures to be verified by DP4. At least one\
    is required, if automatic diastereomer generation is used.")

    parser.add_argument("ExpNMR", help="Experimental NMR description, assigned\
    with the atom numbers from the structure file")

    parser.add_argument("-s", "--solvent", help="Specify solvent to use\
    for dft calculations")
    parser.add_argument("-q", "--queue", help="Specify queue for job submission\
    on ziggy", default=settings.queue)
    parser.add_argument("--TimeLimit", help="Specify job time limit for jobs\
    on ziggy or darwin", type=int)

    parser.add_argument("--nProc", help="Specify number of processor cores\
    to use for Gaussian calculations", type=int, default=1)
    parser.add_argument("--batch", help="Specify max number of jobs per batch",
                        type=int, default=settings.MaxConcurrentJobsZiggy)
    parser.add_argument("--project", help="Specify project for job submission\
    on darwin", default=settings.project)
    parser.add_argument("--ConfLimit", help="Specify maximum number of \
    conformers per structure. If above this, adaptive RMSD pruning will be \
    performed", type=int, default=settings.PerStructConfLimit)

    parser.add_argument("--MaxConfE", help="Specify maximum MMFF energy \
    allowed before conformer is discarded before DFT stage", type=float, \
                        default=settings.MaxCutoffEnergy)

    parser.add_argument("-r", "--rot5", help="Manually generate conformers for\
    5-memebered rings", action="store_true")

    parser.add_argument('--ra', help="Specify ring atoms, for the ring to be\
    rotated, useful for molecules with several 5-membered rings")
    parser.add_argument('-S', '--Stats', help="Specify the stats model and\
    parameters")

    parser.add_argument("--AssumeDFTDone", help="Assume RMSD pruning, DFT setup\
    and DFT calculations have been run already", action="store_true")
    parser.add_argument("--AssumeConverged", help="Assume DFT optimizations have" + \
                                                  " converged and can be used for NMR and or energy calcs",
                        action="store_true")
    parser.add_argument("--UseExistingInputs", help="Use previously generated\
    DFT inputs, avoids long conf pruning and regeneration", action="store_true")
    parser.add_argument("--NoConfPrune", help="Skip RMSD pruning, use all\
    conformers in the energy window", action="store_true")

    parser.add_argument('-c', '--StereoCentres', help="Specify\
    stereocentres for diastereomer generation")
    parser.add_argument("--OptCycles", help="Specify max number of DFT geometry\
    optimization cycles", type=int, default=settings.MaxDFTOptCycles)
    parser.add_argument("--OptStep", help="Specify the max step size\
    Gaussian should take in optimization, default is 30", type=int, default=settings.OptStepSize)
    parser.add_argument("--FC", help="Calculate force constants before optimization", action="store_true")

    parser.add_argument('-n', '--Charge', help="Specify\
    charge of the molecule. Do not use when input files have different charges")
    parser.add_argument('-B', '--nBasisSet', help="Selects the basis set for\
    DFT NMR calculations", default=settings.nBasisSet)
    parser.add_argument('-F', '--nFunctional', help="Selects the functional for\
    DFT NMR calculations", default=settings.nFunctional)
    parser.add_argument('--eBasisSet', help="Selects the basis set for\
    DFT energy calculations", default=settings.eBasisSet)
    parser.add_argument('--eFunctional', help="Selects the functional for\
    DFT energy calculations", default=settings.eFunctional)
    parser.add_argument('-f', '--ff', help="Selects force field for the \
    conformational search, implemented options 'mmff' and 'opls' (2005\
    version)", choices=['mmff', 'opls'], default=settings.ForceField)

    parser.add_argument('--OutputFolder', help="Directory for dp4 output default is cwd", default=settings.OutputFolder)

    parser.add_argument('--Smiles', help="txt file input containing smiles strings on separate lines",
                        default=settings.Smiles)

    parser.add_argument('--Smarts', help="txt file input containing smarts strings on separate lines",
                        default=settings.Smarts)

    parser.add_argument('--InChIs', help="txt file input containing inchi strings on separate lines",
                        default=settings.InChIs)

    args = parser.parse_args()
    print(args.StructureFiles)
    print(args.ExpNMR)

    settings.Title = args.ExpNMR
    settings.NMRsource = args.ExpNMR
    settings.Workflow = args.workflow

    settings.DFT = args.dft
    settings.queue = args.queue
    settings.ScriptDir = getScriptPath()
    settings.ForceField = args.ff
    settings.PerStructConfLimit = args.ConfLimit
    settings.MaxCutoffEnergy = args.MaxConfE
    settings.nBasisSet = args.nBasisSet
    settings.nFunctional = args.nFunctional
    settings.eBasisSet = args.eBasisSet
    settings.eFunctional = args.eFunctional
    settings.nProc = args.nProc
    settings.MaxConcurrentJobs = args.batch
    settings.project = args.project
    settings.MaxDFTOptCycles = args.OptCycles
    settings.OptStepSize = args.OptStep
    if args.FC:
        settings.CalcFC = True

    if args.TimeLimit:
        settings.TimeLimit = args.TimeLimit

    if args.Stats is not None:
        settings.StatsModel = (args.Stats)[0]
        settings.StatsParamFile = (args.Stats)[1:]

    settings.MM = args.mm

    if args.StepCount is not None:
        settings.MMstepcount = int(args.StepCount)
    if args.Charge is not None:
        settings.charge = int(args.Charge)
    if args.StereoCentres is not None:
        settings.SelectedStereocentres = \
            [int(x) for x in (args.StereoCentres).split(',')]
    if args.NoConfPrune:
        settings.ConfPrune = False
    if args.AssumeDFTDone:
        settings.AssumeDone = True
    if args.AssumeConverged:
        settings.AssumeConverged = True
    if args.UseExistingInputs:
        settings.UseExistingInputs = True
    if args.solvent:
        settings.Solvent = args.solvent
    if args.rot5:
        settings.Rot5Cycle = True
    if args.ra is not None:
        settings.RingAtoms = \
            [int(x) for x in (args.ra).split(',')]

    if settings.StatsParamFile != 'none':
        if os.path.isfile(settings.StatsParamFile):
            print("Statistical parameter file found at " + settings.StatsParamFile)
        elif (not os.path.isfile(settings.StatsParamFile)) and \
                os.path.isfile(settings.ScriptDir + settings.StatsParamFile):
            settings.StatsParamFile = settings.ScriptDir + settings.StatsParamFile
            print("Statistical parameter file found at " + settings.StatsParamFile)
        elif (not os.path.isfile(settings.StatsParamFile)) and \
                (not os.path.isfile(settings.ScriptDir + settings.StatsParamFile)):
            print("Stats file not found, quitting.")

    now = datetime.datetime.now()
    settings.StartTime = now.strftime('%d%b%H%M')

    settings.user = getpass.getuser()
    settings.DarwinScrDir.replace('/u/', settings.user)

    with open('cmd.log', 'a') as f:
        f.write(' '.join(sys.argv) + '\n')

    settings.InputFiles = args.StructureFiles

    settings.Smiles = args.Smiles
    settings.Smarts = args.Smarts
    settings.InChIs = args.InChIs

    settings.NMRsource = args.ExpNMR

    NMR_files(settings.NMRsource)

    # check if NMR data has been passed from the cwd or the full path

    settings.OutputFolder = Path(args.OutputFolder)

    main(settings)


if __name__ == '__main__':
    run()

