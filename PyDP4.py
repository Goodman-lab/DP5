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

import NMR
import Tinker
import MacroModel
import DP4n as DP4
import sys
import os
import datetime
import argparse
import importlib

DFTpackages = [['n', 'w', 'j', 'g', 'z', 'd'],['NWChem', 'NWChemZiggy', 'Jaguar', 'Gaussian', 'GaussianZiggy', 'GaussianDarwin']]
'''
if os.name == 'nt':
    import pyximport
    pyximport.install()
    import ConfPrune
else:
    import pyximport
    pyximport.install()
    import ConfPrune
'''
class Paths:
    TinkerPath = '~/tinker7/bin/scan '
    SCHRODINGER = '/usr/local/shared/schrodinger/current'

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
    AssumeConverged = False     # Assume all optimizations have converged, do NMR and/or energy calcs on existing DFT geometries
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
    ConfPrune = False           # Should we prune conformations?
    PerStructConfLimit = 100    # Max numbers of conformers allowed per structure for DFT stages
    InitialRMSDcutoff = 0.75    # Initial RMSD threshold for pruning
    MaxCutoffEnergy = 10.0      # Max conformer MM energy in kJ/mol to allow

    # --- DFT ---
    MaxDFTOptCycles = 50        # Max number of DFT geometry optimization cycles to request.
    CalcFC = False              # Calculate QM force constants before optimization
    OptStepSize = 30            # Max step Gaussian should take in geometry optimization
    charge = None               # Manually specify charge for DFT calcs
    nBasisSet = "6-311g(d)"     # Basis set for NMR calcs
    nFunctional = "mPW1PW91"    # Functional for NMR calcs
    oBasisSet = "6-31g(d,p)"    # Basis set for geometry optimizations
    oFunctional = "b3lyp"       # Functional for geometry optimizations
    eBasisSet = "def2tzvp"      # Basis set for energy calculations
    eFunctional = "m062x"       # Functional for energy calculations

    # --- Computational clusters ---
    """ These should probably be moved to relevant *.py files as Cambridge specific """
    user = 'ah809'              # Linux user on computational clusters, not used for local calcs
    TimeLimit = 24              # Queue time limit on comp clusters
    queue = 'SWAN'              # Which queue to use on Ziggy
    project = 'GOODMAN-SL3-CPU' # Which project to use on Darwin
    DarwinScrDir = '/home/ah809/rds/hpc-work/'  # Which scratch directory to use on Darwin
    StartTime = ''              # Automatically set on launch, used for folder names
    nProc = 1                   # Cores used per job, must be less than node size on cluster
    DarwinNodeSize = 32         # Node size on current CSD3
    MaxConcurrentJobsZiggy = 75      # Max concurrent jobs to submit on ziggy
    MaxConcurrentJobsDarwin = 320 # Max concurrent jobs to submit on CSD3

    # --- NMR analysis ---
    TMS_SC_C13 = 191.69255      # Default TMS reference C shielding constant (from B3LYP/6-31g**)
    TMS_SC_H1 = 31.7518583      # Default TMS reference H shielding constant (from B3LYP/6-31g**)

    # --- Stats ---
    StatsModel = 'g'            # What statistical model type to use
    StatsParamFile = ''         # Where to find statistical model parameters

settings = Settings()

# Data structure keeping all of isomer data in one place.
class Isomer:
    def __init__(self, InputFile, Charge=-100):
        self.InputFile = InputFile  # Initial structure input file
        self.BaseName = InputFile # Basename for other files
        self.Atoms = []             # Element labels
        self.Conformers = []        # from conformational search, list of atom coordinate lists
        self.MMCharge = 0           # charge from conformational search
        self.ExtCharge = Charge     # externally provided charge
        self.RMSDCutoff = 0         # RMSD cutoff eventually used to get the conformer number below the limit
        self.DFTConformers = []     # from DFT optimizations, list of atom coordinate lists
        self.ConfIndices = []       # List of conformer indices from the original conformational search for reference
        self.MMEnergies = []        # Corresponding MM energies in kj/mol
        self.DFTEnergies = []       # Corresponding DFT energies in hartrees
        self.Energies = []          # Final energies used in conformer population prediction in kj/mol
        self.Populations = []       # Conformer populations
        self.OptInputFiles = []     # list of DFT NMR input file names
        self.OptOutputFiles = []    # list of DFT NMR output file names
        self.EInputFiles = []     # list of DFT NMR input file names
        self.EOutputFiles = []    # list of DFT NMR output file names
        self.NMRInputFiles = []     # list of DFT NMR input file names
        self.NMROutputFiles = []    # list of DFT NMR output file names
        self.ShieldingLabels = []   # A list of atom labels corresponding to the shielding values
        self.ConformerShieldings = [] # list of calculated NMR shielding constant lists for every conformer
        self.BoltzmannShieldings = []  # Boltzmann weighted NMR shielding constant list for the isomer
        self.Cshifts = []           # Calculated C NMR shifts
        self.Hshifts = []           # Calculated H NMR
        self.Clabels = []
        self.Hlabels = []
        self.Cexp = []              # Experimental C NMR shifts
        self.Hexp = []              # Experimental H NMR shifts

def main(settings):

    print("==========================")
    print("PyDP4 script,\nintegrating Tinker/MacroModel,")
    print("Gaussian/NWChem/Jaguar and DP4\nv1.0")
    print("\nCopyright (c) 2015-2019 Kristaps Ermanis, Alexander Howarth, Jonathan M. Goodman")
    print("Distributed under MIT license")
    print("==========================\n\n")

    print("Initial input files: " + str(settings.InputFiles))
    print("NMR file: " + str(settings.NMRsource))
    print("Workflow: " + str(settings.Workflow))

    # Check the number of input files, generate some if necessary
    if ('g' in settings.Workflow) and len(settings.InputFiles) == 1:
        import InchiGen
        print("\nGenerating diastereomers...")
        settings.InputFiles = InchiGen.GenDiastereomers(settings.InputFiles[0], settings.SelectedStereocentres)

    print("Generated input files: " + str(settings.InputFiles) + '\n')

    # Create isomer data structures
    Isomers = [Isomer(f) for f in settings.InputFiles]

    # Run conformational search, if requested
    if ('m' in settings.Workflow) and not(settings.AssumeDone or settings.UseExistingInputs):
        if settings.MM == 't':
            print('\nSetting up Tinker files...')
            print('Tinker support not implemented yet!')
            quit()

            TinkerInputs = Tinker.SetupTinker(settings)

            print('\nRunning Tinker...')
            TinkerOutputs = Tinker.RunTinker(TinkerInputs, settings)

            Isomers = Tinker.ReadConformers(TinkerOutputs, Isomers)

        elif settings.MM == 'm':
            print('\nSetting up MacroModel files...')
            MacroModelInputs = MacroModel.SetupMacroModel(settings)
            print("MacroModel inputs: " + str(MacroModelInputs))
            print('Running MacroModel...')
            MacroModelOutputs = MacroModel.RunMacroModel(MacroModelInputs, settings)
            print('\nReading conformers...')
            Isomers = MacroModel.ReadConformers(MacroModelOutputs, Isomers, settings)
            print('Energy window: ' + str(settings.MaxCutoffEnergy) + ' kJ/mol')
            for iso in Isomers:
                print(iso.InputFile + ": "+ str(len(iso.Conformers)) + ' conformers read within energy window' )
    else:
        print('No conformational search was requested. Skipping...')
        settings.ConfPrune = False

    # Prune conformations, if requested.
    # For each isomer, the conformers list is replaced with a smaller list of conformers
    if (settings.ConfPrune) and not(settings.AssumeDone or settings.UseExistingInputs):
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

    if not(settings.AssumeDone):

        # Run DFT optimizations, if requested
        if ('o' in settings.Workflow):

            now = datetime.datetime.now()
            settings.StartTime = now.strftime('%d%b%H%M')

            print('\nSetting up geometry optimization calculations...')
            Isomers = DFT.SetupOptCalcs(Isomers, settings)
            print('\nRunning geometry optimization calculations...')
            Isomers = DFT.RunOptCalcs(Isomers, settings)

            print('\nReading DFT optimized geometries...')

            Isomers = DFT.ReadGeometries(Isomers)

            #Add convergence check here before continuing with calcs!
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
            Isomers = DFT.ReadDFTGeometries(Isomers, settings)

        # Read DFT single-point energies, if requested
        if ('e' in settings.Workflow):
            Isomers = DFT.ReadDFTEnergies(Isomers, settings)

        # Read DFT NMR data, if requested
        if ('n' in settings.Workflow):
            Isomers = DFT.ReadShieldings(Isomers)
            Isomers = DFT.ReadDFTEnergies(Isomers)


    if not(NMR.NMRDataValid(Isomers)) or ('n' not in settings.Workflow):

        print('\nNo NMR data calculated, quitting...')
        quit()

    if 's' in settings.Workflow:

        print('\nSetting TMS computational NMR shielding constant references')
        settings.TMS_SC_C13, settings.TMS_SC_H1 = NMR.GetTMSConstants(settings)

        print('\nConverting DFT data to NMR shifts...')
        Isomers = NMR.CalcBoltzmannWeightedShieldings(Isomers)
        Isomers = NMR.CalcNMRShifts(Isomers, settings)

        print('\nReading experimental NMR data...')
        NMRData = NMR.NMRData(settings)

        '''
        print("Conformation data:")
        NMR.PrintConformationData(Isomers)
        '''

        if NMRData.Type == 'desc':

            print('Experimental NMR description found and read.')

            # performs a pairwise assignment

            Isomers = NMR.PairwiseAssignment(Isomers,NMRData)

            # print('Clabels: ' + str(NMRData.Clabels))
            print('Cshifts: ' + str(NMRData.Cshifts))

            # print('Hlabels: ' + str(NMRData.Hlabels))
            print('Hshifts: ' + str(NMRData.Hshifts))

            print('Equivalents: ' + str(NMRData.Equivalents))
            print('Omits: ' + str(NMRData.Omits))

        elif NMRData.Type == 'fid':

            if os.path.exists(str(settings.NMRsource) + "/Proton"):

                from Proton_assignment import AssignProton

                from Proton_plotting import PlotProton

                Isomers = AssignProton(NMRData,Isomers,settings)

                PlotProton(NMRData, Isomers, settings)

            if os.path.exists(str(settings.NMRsource) + "/Carbon"):

                from Carbon_assignment import AssignCarbon

                from Carbon_plotting import PlotCarbon

                Isomers = AssignCarbon(NMRData,Isomers,settings)

                PlotCarbon(NMRData, Isomers, settings)

            print('Raw FID NMR datafound and read.')


        #print('\nProcessing experimental NMR data...')

        #NMRdata = NMR.ProcessNMRData(Isomers, settings.NMRsource, settings)

        print('\nCalculating DP4 probabilities...')

        DP4data = DP4.DP4data()

        DP4data = DP4.ProcessIsomers(DP4data,Isomers)

        DP4data = DP4.InternalScaling(DP4data)

        DP4data = DP4.CalcProbs(DP4data,settings)

        DP4data = DP4.CalcDP4(DP4data)

        DP4data = DP4.MakeOutput(DP4data,Isomers,settings)

        quit()

        #print(DP4.FormatData(DP4data))


    else:
        print('\nNo DP4 analysis requested.')

    print('\nPyDP4 process completed successfully.')


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


if __name__ == '__main__':

    #Read config file and fill in settings in from that
    #These are then overridden by any explicit parameters given through the command line

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
    parser.add_argument("--AssumeConverged", help="Assume DFT optimizations have" + \
    " converged and can be used for NMR and or energy calcs", action="store_true")
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
        settings.SelectedStereocentres =\
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

    now = datetime.datetime.now()
    settings.StartTime = now.strftime('%d%b%H%M')

    with open('cmd.log', 'a') as f:
        f.write(' '.join(sys.argv) + '\n')

    settings.InputFiles = args.StructureFiles

    settings.NMRsource = args.ExpNMR
    
    main(settings)
