#!/usr/bin/env python3
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

import NMRAnalysis
import Tinker
import MacroModel

import sys
import os
import datetime
import argparse
import math
import importlib

DFTpackages = [['n', 'w', 'j', 'g', 'z', 'd'],
    ['NWChem.py', 'NWChemZiggy.py', 'Jaguar.py', 'Gaussian.py', 'GaussianZiggy.py', 'GaussianDarwin.py']]

if os.name == 'nt':
    import pyximport
    pyximport.install()
    import ConfPrune
else:
    import pyximport
    pyximport.install()
    import ConfPrune

class Paths:
    TinkerPath = '~/tinker7/bin/scan '
    OBPath = '/home/ke291/Tools/openbabel-install/lib/python2.7/site-packages/'
    SCHRODINGER = '/usr/local/shared/schrodinger/current'
    MOPAC = '/home/ke291/MOPAC/MOPAC2016.exe'

#Assigning the config default values
#Settings are defined roughly in the order they are used in the script
class Settings:
    # --- Main options ---
    MM = 'm'        # m for MacroModel, t for Tinker
    DFT = 'z'       # n, j, g, z or for NWChem, Jaguar or Gaussian
    Workflow = 'gmns' # defines which steps to include in the workflow
                    # g for generate diastereomers
                    # m for molecular mechanics conformational search
                    # o for DFT optimization
                    # e for DFT single-point energies
                    # n for DFT NMR calculation
                    # s for computational and experimental NMR data extraction and stats analysis
    Solvent = ''    # solvent for DFT optimization and NMR calculation
    ScriptDir = ''  # Script directory, automatically set on launch
    InputFiles = [] # Structure input files - can be MacroModel *-out.mae or *sdf files
    NMRsource = ''  # File or folder containing NMR description or data
    Title = 'DP4molecule'       # Title of the calculation, set to NMR file name by default on launch
    AssumeDone = False          # Assume all computations done, only read DFT output data and analyze (use for reruns)
    UseExistingInputs = False   # Don't regenerate DFT inputs, use existing ones. Good for restarting a failed calc

    # --- Diastereomer generation ---
    SelectedStereocentres = []  # which stereocentres to vary for diastereomer generation

    # --- Molecular mechanics ---
    ForceField = 'mmff'         # ff to use for conformational search
    MMstepcount = 10000         # Max number of MM steps to do, if less than MMfactor*rotable_bonds
    MMfactor = 2500             # MMfactor*rotable_bonds gives number of steps to do if less than MMstepcount
    Rot5Cycle = False           # Special dealing with 5-membered saturated rings, see FiveConf.py
    RingAtoms = []              # Define the 5-membered ring, useful if several are present in molecule
    SCHRODINGER = ''            # Define the root folder for Schrodinger software

    # --- Conformer pruning ---
    HardConfLimit = 10000       # Immediately stop if conformers to run exceed this number
    ConfPrune = True            # Should we prune conformations?
    PerStructConfLimit = 100    # Max numbers of conformers allowed per structure for DFT stages
    InitialRMSDcutoff = 0.75    # Initial RMSD threshold for pruning
    MaxCutoffEnergy = 10.0      # Max conformer MM energy in kJ/mol to allow

    # --- DFT ---
    MaxDFTOptCycles = 50        # Max number of DFT geometry optimization cycles to request.
    charge = None               # Manually specify charge for DFT calcs
    nBasisSet = "6-311g(d)"     # Basis set for NMR calcs
    nFunctional = "mPW1PW91"    # Functional for NMR calcs
    oBasisSet = "6-31g(d,p)"    # Basis set for geometry optimizations
    oFunctional = "b3lyp"       # Functional for geometry optimizations
    eBasisSet = "def2tzvp"      # Basis set for energy calculations
    eFunctional = "m062x"       # Functional for energy calculations

    # --- Computational clusters ---
    """ These should probably be moved to relevant *.py files as Cambridge specific """
    user = 'ke291'              # Linux user on computational clusters, not used for local calcs
    TimeLimit = 24              # Queue time limit on comp clusters
    queue = 's1'                # Which queue to use
    DarwinScrDir = '/home/ke291/rds/hpc-work/'  # Which scratch directory to use on Darwin
    StartTime = ''              # Automatically set on launch, used for folder names
    nProc = 1                   # Cores used per job, must be less than node size on cluster
    DarwinNodeSize = 32         # Node size on current CSD3
    MaxConcurrentJobs = 75      # Max concurrent jobs to submit on ziggy
    MaxConcurrentJobsDarwin = 320 # Max concurrent jobs to submit on CSD3

    # --- NMR analysis ---
    TMS_SC_C13 = 191.69255      #Default TMS reference C shielding constant (from B3LYP/6-31g**)
    TMS_SC_H1 = 31.7518583      #Default TMS reference H shielding constant (from B3LYP/6-31g**)

    # --- Stats ---
    StatsModel = 'g'            #What statistical model type to use
    StatsParamFile = ''         #Where to find statistical model parameters

settings = Settings()

# Data struture keeping all of isomer data in one place.
class Isomer:
    def __init__(self, InputFile):
        self.InputFile = InputFile  # Initial structure input file
        self.atoms = []             # Element labels
        self.Conformers = []        # from conformational search, list of atom coordinate lists
        self.DFTConformers = []     # from DFT optimizations, list of atom coordinate lists
        self.ConfIndices = []       # List of conformer indices from the original conformational search for reference
        self.MMEnergies = []        # Corresponding MM energies
        self.DFTEnergies = []       # Corresponding DFT energies
        self.OptInputFiles = []     # list of DFT NMR input file names
        self.OptOutputFiles = []    # list of DFT NMR output file names
        self.EInputFiles = []     # list of DFT NMR input file names
        self.EOutputFiles = []    # list of DFT NMR output file names
        self.NMRInputFiles = []     # list of DFT NMR input file names
        self.NMROutputFiles = []    # list of DFT NMR output file names
        self.Cshifts = []
        self.Hshifts = []


def main(settings):

    print("==========================")
    print("PyDP4 script,\nintegrating Tinker/MacroModel,")
    print("Gaussian/NWChem/Jaguar and DP4\nv1.0")
    print("\nCopyright (c) 2015-2019 Kristaps Ermanis, Alexander Howarth, Jonathan M. Goodman")
    print("Distributed under MIT license")
    print("==========================\n\n")

    print(settings.InputFiles)
    # Check the number of input files, generate some if necessary
    if ('g' in settings.Workflow) and len(settings.InputFiles) == 1:
        import InchiGen
        settings.InputFiles = InchiGen.GenDiastereomers(settings.InputFiles[0], settings.SelectedStereocentres)

    print("Input files: " + str(settings.InputFiles))
    print("NMR file: " + str(settings.NMRsource))

    # Create isomer data structures
    Isomers = [Isomer(f) for f in settings.InputFiles]

    # Run conformational search, if requested
    if ('m' in settings.Workflow) and not(settings.AssumeDone or settings.UseExistingInputs):
        if settings.MM == 't':
            print('\nSetting up Tinker files...')
            TinkerInputs = Tinker.SetupTinker(settings)

            print('\nRunning Tinker...')
            TinkerOutputs = Tinker.RunTinker(TinkerInputs, settings)

            Isomers = Tinker.ReadConformers(TinkerOutputs, Isomers)

        elif settings.MM == 'm':
            print('\nSetting up MacroModel files...')
            MacroModelInputs = MacroModel.SetupMacroModel(settings)

            print('\nRunning Tinker...')
            MacroModelOutputs = MacroModel.RunMacroModel(MacroModelInputs, settings)

            Isomers = MacroModel.ReadConformers(MacroModelOutputs, Isomers)

    else:
        print('No conformational search was requested. Skipping...')
        settings.ConfPrune = False

    # Prune conformations, if requested.
    # For each isomer, the conformers list is replaced with a smaller list of conformers
    if (settings.ConfPrune) and not(settings.AssumeDone or settings.UseExistingInputs):
        Isomers = ConfPrune.AdaptRMSDPrune(Isomers, settings)

    if ('d' in settings.Workflow) or ('o' in settings.Workflow) \
            or ('e' in settings.Workflow) or settings.AssumeDone:
        DFT = ImportDFT(settings.DFT)

    if not(settings.AssumeDone):

        # Run DFT optimizations, if requested
        if ('o' in settings.Workflow):
            Isomers = DFT.SetupOptCalcs(Isomers, settings)
            DFT.RunOptCalcs(Isomers, settings)
            Isomers = DFT.ReadDFTGeometries(Isomers, settings)

        # Run DFT single-point energy calculations, if requested
        if ('e' in settings.Workflow):
            Isomers = DFT.SetupECalcs(Isomers, settings)
            DFT.RunECalcs(Isomers, settings)
            Isomers = DFT.ReadDFTEnergies(Isomers, settings)

        # Run DFT NMR calculations, if requested
        if ('d' in settings.Workflow):
            Isomers = DFT.SetupNMRCalcs(Isomers, settings)
            DFT.RunNMRCalcs(Isomers, settings)
            SetTMSConstants(settings)
            Isomers = DFT.ReadNMRShifts(Isomers, settings)

    else:
        # Run DFT optimizations, if requested
        if ('o' in settings.Workflow):
            Isomers = DFT.ReadDFTGeometries(Isomers, settings)

        # Run DFT single-point energy calculations, if requested
        if ('e' in settings.Workflow):
            Isomers = DFT.ReadDFTEnergies(Isomers, settings)

        # Run DFT NMR calculations, if requested
        if ('d' in settings.Workflow):
            SetTMSConstants(settings)
            Isomers = DFT.ReadNMRShifts(Isomers, settings)


    if not(NMRDataValid(Isomers)) or ('n' not in settings.Workflow):

        print('No NMR data calculated, quitting...')
        quit()

    NMRAnalysis.main(Isomers, settings)

    print('DP4 process completed successfully.')


def NMRDataValid(Isomers):

    for isomer in Isomers:
        if (len(isomer.Cshifts) == 0) and (len(isomer.Hshifts) == 0):
            return False

    return True


# Selects which DFT package to import, returns imported module
def ImportDFT(dft):
    if dft in DFTpackages[0]:
        DFTindex = DFTpackages[0].index(dft)
        DFT = importlib.import_module(DFTpackages[1][DFTindex])
    else:
        print("Invalid DFT package selected")
        quit()

    return DFT


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def SetTMSConstants(settings):
    
    TMSfile = open(settings.ScriptDir + '/TMSdata', 'r')
    TMSdata = TMSfile.readlines()
    TMSfile.close()
    
    for i, line in enumerate(TMSdata):
        buf = line.split(' ')
        if len(buf)>1:
            if settings.Solvent != '':
                if buf[0].lower() == settings.Functional.lower() and \
                   buf[1].lower() == settings.BasisSet.lower() and \
                   buf[2].lower() == settings.Solvent.lower():
                    
                    print("Setting TMS references to " + buf[3] + " and " + \
                        buf[4] + "\n")
                    settings.TMS_SC_C13 = float(buf[3])
                    settings.TMS_SC_H1 = float(buf[4])
                    return
            else:
                if buf[0].lower() == settings.Functional.lower() and \
                   buf[1].lower() == settings.BasisSet.lower() and \
                   buf[2].lower() == 'none':
                    
                    print("Setting TMS references to " + buf[3] + " and " + \
                        buf[4] + "\n")
                    settings.TMS_SC_C13 = float(buf[3])
                    settings.TMS_SC_H1 = float(buf[4])
                    return
    
    print("No TMS reference data found for these conditions, using defaults\n")
    print("Unscaled shifts might be inaccurate, use of unscaled models is" + \
        " not recommended.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PyDP4 script to setup\
    and run Tinker, Gaussian (on ziggy) and DP4')
    parser.add_argument('-w', '--workflow', help="Defines which steps to include in the workflow, " +
    "can contain g for generate diastereomers, m for molecular mechanics conformational search, " +
    "o for DFT optimization, e for DFT single-point energies, n for DFT NMR calculation, " +
    "s for computational and experimental NMR data extraction and stats analysis, default is 'gmns'", default=settings.Workflow)

    parser.add_argument('-m', '--mm', help="Select molecular mechanics program,\
    t for tinker or m for macromodel, default is m", choices=['t', 'm'],
    default='m')
    parser.add_argument('-d', '--dft', help="Select DFT program, j for Jaguar,\
    g for Gaussian, n for NWChem, z for Gaussian on ziggy, d for Gaussian on \
    Darwin, default is g", choices=DFTpackages[0], default='g')
    parser.add_argument('--StepCount', help="Specify\
    stereocentres for diastereomer generation")
    parser.add_argument('StructureFiles', nargs='+', default=['-'], help=
    "One or more SDF file for the structures to be verified by DP4. At least one\
    is required, if automatic diastereomer generation is used.")
    parser.add_argument("ExpNMR", help="Experimental NMR description, assigned\
    with the atom numbers from the structure file")
    parser.add_argument("-s", "--solvent", help="Specify solvent to use\
    for dft calculations")
    parser.add_argument("-q", "--queue", help="Specify queue for job submission\
    on ziggy", default='s1')
    parser.add_argument("--TimeLimit", help="Specify job time limit for jobs\
    on ziggy or darwin", type=int)

    parser.add_argument("--nProc", help="Specify number of processor cores\
    to use for Gaussian calculations", type=int, default=1)
    parser.add_argument("--batch", help="Specify max number of jobs per batch",
    type=int, default=settings.MaxConcurrentJobs)
    parser.add_argument("--ConfLimit", help="Specify maximum number of \
    conformers per structure. If above this, adaptive RMSD pruning will be \
    performed", type=int, default=settings.PerStructConfLimit)

    parser.add_argument("--MaxConfE", help="Specify maximum MMFF energy \
    allowed before conformer is discarded before DFT stage", type=float,\
    default=settings.MaxCutoffEnergy)

    parser.add_argument("-r", "--rot5", help="Manually generate conformers for\
    5-memebered rings", action="store_true")

    parser.add_argument('--ra', help="Specify ring atoms, for the ring to be\
    rotated, useful for molecules with several 5-membered rings")
    parser.add_argument('-S', '--Stats', help="Specify the stats model and\
    parameters")

    parser.add_argument("--AssumeDFTDone", help="Assume RMSD pruning, DFT setup\
    and DFT calculations have been run already", action="store_true")
    parser.add_argument("--UseExistingInputs", help="Use previously generated\
    DFT inputs, avoids long conf pruning and regeneration", action="store_true")
    parser.add_argument("--NoConfPrune", help="Skip RMSD pruning, use all\
    conformers in the energy window", action="store_true")

    parser.add_argument('-c', '--StereoCentres', help="Specify\
    stereocentres for diastereomer generation")
    parser.add_argument('-o', '--DFTOpt', help="Optimize geometries at DFT\
    level before NMR prediction", action="store_true")
    parser.add_argument("--OptCycles", help="Specify max number of DFT geometry\
    optimization cycles", type=int, default=settings.MaxDFTOptCycles)
    parser.add_argument('-n', '--Charge', help="Specify\
    charge of the molecule. Do not use when input files have different charges")
    parser.add_argument('-B', '--BasisSet', help="Selects the basis set for\
    DFT calculations", default=settings.nBasisSet)
    parser.add_argument('-F', '--Functional', help="Selects the functional for\
    DFT calculations", default=settings.nFunctional)
    parser.add_argument('-f', '--ff', help="Selects force field for the \
    conformational search, implemented options 'mmff' and 'opls' (2005\
    version)", choices=['mmff', 'opls'], default=settings.ForceField)
    args = parser.parse_args()
    print(args.StructureFiles)
    print(args.ExpNMR)
    settings.Title = args.ExpNMR
    settings.NMRsource = args.ExpNMR

    settings.DFT = args.dft
    settings.queue = args.queue
    settings.ScriptDir = getScriptPath()
    settings.ForceField = args.ff
    settings.PerStructConfLimit = args.ConfLimit
    settings.MaxCutoffEnergy = args.MaxConfE
    settings.nBasisSet = args.BasisSet
    settings.nFunctional = args.Functional
    settings.nProc = args.nProc
    settings.MaxConcurrentJobs = args.batch
    settings.MaxDFTOptCycles = args.OptCycles
    
    if args.TimeLimit:
        settings.TimeLimit = args.TimeLimit

    if args.Stats is not None:
        settings.StatsModel = (args.Stats)[0]
        settings.StatsParamFile = (args.Stats)[1:]

    settings.MM = args.mm

    if args.DFTOpt:
        settings.Workflow = settings.Workflow + 'o'
    if args.StepCount is not None:
        settings.MMstepcount = int(args.StepCount)
    if args.Charge is not None:
        settings.charge = int(args.Charge)
    if args.StereoCentres is not None:
        settings.SelectedStereocentres =\
            [int(x) for x in (args.StereoCentres).split(',')]
    if args.NoConfPrune:
        settings.ConfPrune = False
    if args.AssumeDFTDone:
        settings.AssumeDone = True
    if args.UseExistingInputs:
        settings.UseExistingInputs = True
    if args.solvent:
        settings.Solvent = args.solvent
    if args.rot5:
        settings.Rot5Cycle = True
    if args.ra is not None:
        settings.RingAtoms =\
            [int(x) for x in (args.ra).split(',')]
    
    if settings.StatsParamFile != '':
        if os.path.isfile(settings.StatsParamFile):
            print("Statistical parameter file found at " + settings.StatsParamFile)
        elif (not os.path.isfile(settings.StatsParamFile)) and\
            os.path.isfile(settings.ScriptDir+settings.StatsParamFile):
                settings.StatsParamFile = settings.ScriptDir+settings.StatsParamFile
                print("Statistical parameter file found at " + settings.StatsParamFile)
        elif (not os.path.isfile(settings.StatsParamFile)) and\
            (not os.path.isfile(settings.ScriptDir+settings.StatsParamFile)):
            print("Stats file not found, quitting.")
            quit()
    
    now = datetime.datetime.now()
    settings.StartTime = now.strftime('%d%b%H%M')

    with open('cmd.log', 'a') as f:
        f.write(' '.join(sys.argv) + '\n')

    settings.InputFiles = args.StructureFiles

    settings.NMRsource = args.ExpNMR
    
    main(settings)
