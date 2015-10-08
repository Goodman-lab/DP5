#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyDP4 integrated workflow for the running of MM, DFT GIAO calculations and
DP4 analysis
v0.4

Copyright (c) 2015 Kristaps Ermanis, Jonathan M. Goodman

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
Updated on Mon Jul 30 2015

@author: ke291

The main file, that should be called to start the PyDP4 workflow.
Interprets the arguments and takes care of the general workflow logic.
"""

from __future__ import division

import Gaussian
import NMRAnalysis
import Tinker
import MacroModel
import NWChem

import glob
import sys
import os
import datetime
import argparse
import math


#Assigning the config default values
class Settings:
    MMTinker = True
    MMMacromodel = False
    DFT = 'z'
    Rot5Cycle = False
    RingAtoms = []
    GenDS = True
    GenTaut = False
    GenProt = False
    Solvent = ''
    DFTOpt = False
    PDP4 = True
    EP5 = False
    jKarplus = False
    jFC = False
    jJ = False
    queue = 's1'
    TinkerPath = '~/tinker7/bin/scan '
    OBPath = '/home/ke291/Tools/openbabel-install/lib/python2.7/site-packages/'
    SCHRODINGER = '/usr/local/shared/schrodinger/current'
    ScriptDir = ''
    user = 'ke291'
    MMstepcount = 10000
    MMfactor = 500  # nsteps = MMfactor*degrees of freedom
    HardConfLimit = 10000
    MaxConcurrentJobs = 75
    PerStructConfLimit = 100
    InitialRMSDcutoff = 0.75
    MaxCutoffEnergy = 10.0
    NTaut = 1
    LogFile = 'PyDP4.log'
    AssumeDone = False
    GenOnly = False
    SelectedStereocentres = []
    charge = None
    BasicAtoms = []
    ForceField = 'mmff'
    BasisSet = "6-31g(d,p)"
    Functional = "b3lyp"

settings = Settings()


def main(filename, ExpNMR, nfiles):

    print "=========================="
    print "PyDP4 script,\nintegrating Tinker/MacroModel,"
    print "Gaussian/NWChem and DP4\nv0.4"
    print "\nCopyright (c) 2015 Kristaps Ermanis, Jonathan M. Goodman"
    print "Distributed under MIT license"
    print "==========================\n\n"

    if nfiles < settings.NTaut or nfiles % settings.NTaut != 0:
        print "Invalid number of tautomers/input files - number of input files\
        must be a multiple of number of tautomers"
        quit()

    #Check the number of input files, generate some if necessary
    if nfiles == 1:
        import InchiGen
        if len(settings.SelectedStereocentres) > 0:
            numDS, inpfiles = InchiGen.GenSelectDiastereomers(filename,
                                                settings.SelectedStereocentres)
        else:
            numDS, inpfiles = InchiGen.GenDiastereomers(filename)
        if settings.GenTaut:
            newinpfiles = []
            for ds in inpfiles:
                print "Generating tautomers for " + ds
                settings.NTaut, files = InchiGen.GenTautomers(ds)
                newinpfiles.extend(files)
            inpfiles = list(newinpfiles)
        if settings.GenProt:
            newinpfiles = []
            for ds in inpfiles:
                print "Generating protomers for " + ds
                settings.NTaut, files = InchiGen.GenProtomers(ds,
                                                        settings.BasicAtoms)
                newinpfiles.extend(files)
            inpfiles = list(newinpfiles)
    else:
        inpfiles = filename
        if settings.GenTaut:
            numDS = nfiles
            import InchiGen
            newinpfiles = []
            for ds in inpfiles:
                print "Generating tautomers for " + ds
                settings.NTaut, files = InchiGen.GenTautomers(ds)
                newinpfiles.extend(files)
            inpfiles = list(newinpfiles)
        else:
            numDS = int(nfiles/settings.NTaut)
            if numDS == 1:
                import InchiGen
                for f in filename:
                    tdiastereomers = []
                    numDS, tinpfiles = InchiGen.GenDiastereomers(f)
                    tdiastereomers.append(tinpfiles)
                tinpfiles = zip(*tdiastereomers)
                inpfiles = []
                for ds in tinpfiles:
                    inpfiles.extend(list(ds))

    print inpfiles

    #Check the existence of mm output files
    MMRun = False

    if settings.MMTinker:
        #Check if there already are Tinker output files with the right names
        tinkfiles = glob.glob('*.tout')
        mminpfiles = []
        for f in inpfiles:
            if f + '.tout' in tinkfiles and (f + 'rot.tout' in
                                    tinkfiles or settings.Rot5Cycle is False):
                if len(mminpfiles) == 0:
                    MMRun = True
            else:
                MMRun = False
                mminpfiles.append(f)
    else:
        #Check if there already are MacroModel output files with the right names
        mmfiles = glob.glob('*.log')
        mminpfiles = []
        for f in inpfiles:
            if f + '.log' in mmfiles and (f + 'rot.log' in
                                        mmfiles or settings.Rot5Cycle is False):
                if len(mminpfiles) == 0:
                    MMRun = True
            else:
                MMRun = False
                mminpfiles.append(filename)

    if MMRun:
        print 'Conformation search has already been run for these inputs.\
                \nSkipping...'
        if settings.GenOnly:
            print "Input files generated, quitting..."
            quit()
    else:
        if settings.MMTinker:
            print 'Some Tinker files missing.'
            print '\nSeting up Tinker files...'
            Tinker.SetupTinker(len(inpfiles), settings, *mminpfiles)
            if settings.GenOnly:
                print "Input files generated, quitting..."
                quit()
            print '\nRunning Tinker...'
            Tinker.RunTinker(len(inpfiles), settings, *mminpfiles)
        else:
            print 'Some Macromodel files missing.'
            print '\nSetting up Macromodel files...'
            MacroModel.SetupMacromodel(len(inpfiles), settings, *mminpfiles)
            if settings.GenOnly:
                print "Input files generated, quitting..."
                quit()
            print '\nRunning Macromodel...'
            MacroModel.RunMacromodel(len(inpfiles), settings, *mminpfiles)

    if not settings.AssumeDone:
        if settings.DFT == 'z' or settings.DFT == 'g':
            adjRMSDcutoff = Gaussian.AdaptiveRMSD(inpfiles[0], settings)
        elif settings.DFT == 'n' or settings.DFT == 'w':
            adjRMSDcutoff = NWChem.AdaptiveRMSD(inpfiles[0], settings)
        print 'RMSD cutoff adjusted to ' + str(adjRMSDcutoff)
        #Run NWChem setup script for every diastereomer
        print '\nRunning DFT setup...'
        i = 1
        for ds in inpfiles:
            if settings.DFT == 'z' or settings.DFT == 'g':
                print "\nGaussian setup for file " + ds + " (" + str(i) +\
                    " of " +  str(len(inpfiles)) + ")"
                Gaussian.SetupGaussian(ds, ds + 'ginp', 3, settings,
                                       adjRMSDcutoff)
            elif settings.DFT == 'n' or 'w':
                print "\nNWChem setup for file " + ds +\
                    " (" + str(i) + " of " + str(len(inpfiles)) + ")"
                NWChem.SetupNWChem(ds, ds + 'nwinp', 3, settings,
                                   adjRMSDcutoff)
            i += 1
        QRun = False
    else:
        QRun = True

    if settings.DFT == 'z' or settings.DFT == 'g':
        Files2Run = Gaussian.GetFiles2Run(inpfiles, settings)
    elif settings.DFT == 'n' or 'w':
        Files2Run = NWChem.GetFiles2Run(inpfiles, settings)
    print Files2Run
    if len(Files2Run) == 0:
        QRun = True

    if len(Files2Run) > settings.HardConfLimit:
        print "Hard conformation count limit exceeded, DFT calculations aborted."
        quit()

    if QRun:
        print 'DFT has already been run for these inputs. Skipping...'
    else:
        if settings.DFT == 'z':
            print '\nRunning Gaussian on Ziggy...'

            #Run Gaussian jobs on Ziggy cluster in folder named after date
            #and time in the short 1processor job queue
            #and wait until the last file is completed
            now = datetime.datetime.now()
            MaxCon = settings.MaxConcurrentJobs
            if settings.DFTOpt:
                for i in range(len(Files2Run)):
                    Files2Run[i] = Files2Run[i][:-5] + '.com'
            if len(Files2Run) < MaxCon:
                Gaussian.RunOnZiggy(now.strftime('%d%b%H%M'), settings.queue,
                                    Files2Run, settings)
            else:
                print "The DFT calculations will be done in " +\
                    str(math.ceil(len(Files2Run)/MaxCon)) + " batches"
                i = 0
                while (i+1)*MaxCon < len(Files2Run):
                    print "Starting batch nr " + str(i+1)
                    Gaussian.RunOnZiggy(now.strftime('%d%b%H%M')+str(i+1),
                        settings.queue, Files2Run[(i*MaxCon):((i+1)*MaxCon)], settings)
                    i += 1
                print "Starting batch nr " + str(i+1)
                Gaussian.RunOnZiggy(now.strftime('%d%b%H%M')+str(i+1),
                    settings.queue, Files2Run[(i*MaxCon):], settings)

        elif settings.DFT == 'n':
            print '\nRunning NWChem locally...'
            NWChem.RunNWChem(Files2Run, settings)
        elif settings.DFT == 'w':
            print '\nRunning NWChem on Ziggy...'

            #Run NWChem jobs on Ziggy cluster in folder named after date
            #and time in the short 1 processor job queue
            #and wait until the last file is completed
            now = datetime.datetime.now()
            MaxCon = settings.MaxConcurrentJobs
            if len(Files2Run) < MaxCon:
                NWChem.RunOnZiggy(now.strftime('%d%b%H%M'), settings.queue,
                                  Files2Run, settings)
            else:
                print "The DFT calculations will be done in " +\
                    str(math.ceil(len(Files2Run)/MaxCon)) + " batches"
                i = 0
                while (i+1)*MaxCon < len(Files2Run):
                    print "Starting batch nr " + str(i+1)
                    NWChem.RunOnZiggy(now.strftime('%d%b%H%M')+str(i+1),
                        settings.queue, Files2Run[(i*MaxCon):((i+1)*MaxCon)], settings)
                    i += 1
                print "Starting batch nr " + str(i+1)
                NWChem.RunOnZiggy(now.strftime('%d%b%H%M')+str(i+1),
                    settings.queue, Files2Run[(i*MaxCon):], settings)

    if numDS < 2:
        print "DP4 requires at least 2 candidate structures!"
    else:
        allargs = []
        for i in range(numDS):
            allargs.append(settings.NTaut)
            allargs.extend(inpfiles[i*settings.NTaut:(i+1)*settings.NTaut])
        allargs.append(ExpNMR)
        DP4outp = NMRAnalysis.main(numDS, settings, *allargs)
        print '\nWriting the DP4 output to DP4outp'
        DP4_ofile = open(allargs[-1] + '.dp4', 'w')
        DP4_ofile.write(DP4outp)
        DP4_ofile.close()
        print 'DP4 process completed successfully.'


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PyDP4 script to setup\
    and run Tinker, Gaussian (on ziggy) and DP4')
    parser.add_argument('-m', '--mm', help="Select molecular mechanics program,\
    t for tinker or m for macromodel, default is t", choices=['t', 'm'],
    default='t')
    parser.add_argument('-d', '--dft', help="Select DFT program, j for Jaguar,\
    g for Gaussian, n for NWChem, z for Gaussian on ziggy, w for NWChem on\
    ziggy, default is z", choices=['j', 'g', 'n', 'z', 'w'], default='z')
    parser.add_argument('--StepCount', help="Specify\
    stereocentres for diastereomer generation")
    parser.add_argument('StructureFiles', nargs='+', default=['-'], help=
    "One or more SDF file for the structures to be verified by DP4. At least one\
    is required, if automatic diastereomer and tautomer generation is used.\
    One for each candidate structure, if automatic generation is not used")
    parser.add_argument("ExpNMR", help="Experimental NMR description, assigned\
    with the atom numbers from the structure file")
    parser.add_argument("-s", "--solvent", help="Specify solvent to use\
    for dft calculations")
    parser.add_argument("-q", "--queue", help="Specify queue for job submission\
    on ziggy", default='s1')
    parser.add_argument("-t", "--ntaut", help="Specify number of explicit\
    tautomers per diastereomer given in structure files, must be a multiple\
    of structure files", type=int, default=1)
    parser.add_argument("-l", "--ConfLimit", help="Specify maximum number of \
    conformers per structure. If above this, adaptive RMSD pruning will be \
    performed", type=int, default=100)
    parser.add_argument("-r", "--rot5", help="Manually generate conformers for\
    5-memebered rings", action="store_true")
    parser.add_argument("--jJ", help="Calculate coupling constants at DFT\
    level and use in analysis", action="store_true")
    parser.add_argument("--jFC", help="Calculate Fermi contact term of\
    coupling constants at DFT level level and use in analysis", action="store_true")
    parser.add_argument("--jK", help="Use Karplus-type equation to calculate\
    coupling constants and use in analysis", action="store_true")
    parser.add_argument('--ra', help="Specify ring atoms, for the ring to be\
    rotated, useful for molecules with several 5-membered rings")
    parser.add_argument("--AssumeDFTDone", help="Assume RMSD pruning, DFT setup\
    and DFT calculations have been run already", action="store_true")
    parser.add_argument("-g", "--GenOnly", help="Only generate diastereomers\
    and tinker input files, but don't run any calculations", action="store_true")
    parser.add_argument('-c', '--StereoCentres', help="Specify\
    stereocentres for diastereomer generation")
    parser.add_argument('-T', '--GenTautomers', help="Automatically generate\
    tautomers", action="store_true")
    parser.add_argument('-o', '--DFTOpt', help="Optimize geometries at DFT\
    level before NMR prediction", action="store_true")
    parser.add_argument('--ep5', help="Use EP5", action="store_true")
    parser.add_argument('-n', '--Charge', help="Specify\
    charge of the molecule. Do not use when input files have different charges")
    parser.add_argument('-b', '--BasicAtoms', help="Generate protonated states\
    on the specified atoms and consider as tautomers")
    parser.add_argument('-B', '--BasisSet', help="Selects the basis set for\
    DFT calculations", default='6-31g(d,p)')
    parser.add_argument('-F', '--Functional', help="Selects the functional for\
    DFT calculations", default='b3lyp')
    parser.add_argument('-f', '--ff', help="Selects force field for the \
    conformational search, implemented options 'mmff' and 'opls' (2005\
    version)", choices=['mmff', 'opls'], default='mmff')
    args = parser.parse_args()
    print args.StructureFiles
    print args.ExpNMR
    settings.NTaut = args.ntaut
    settings.DFT = args.dft
    settings.queue = args.queue
    settings.ScriptDir = getScriptPath()
    settings.ForceField = args.ff
    settings.PerStructConfLimit = args.ConfLimit
    settings.BasisSet = args.BasisSet
    settings.Functional = args.Functional
    
    if args.jJ:
        settings.jJ = True
        settings.jFC = False
        settings.jKarplus = False
    if args.jFC:
        settings.jFC = True
        settings.jJ = False
        settings.jKarplus = False
    if args.jK:
        settings.jKarplus = True
        settings.jJ = False
        settings.jFC = False
    if args.ep5:
        settings.EP5 = True
        settings.PDP4 = False
    else:
        settings.EP5 = False
        settings.PDP4 = True

    if args.mm == 't':
        settings.MMTinker = True
        settings.MMMacromodel = False
    else:
        settings.MMMacromodel = True
        settings.MMTinker = False
    if args.DFTOpt:
        settings.DFTOpt = True
    if args.BasicAtoms is not None:
        settings.BasicAtoms =\
            [int(x) for x in (args.BasicAtoms).split(',')]
        settings.GenProt = True
    if args.StepCount is not None:
        settings.MMstepcount = int(args.StepCount)
    if args.Charge is not None:
        settings.charge = int(args.Charge)
    if args.GenTautomers:
        settings.GenTaut = True
    if args.StereoCentres is not None:
        settings.SelectedStereocentres =\
            [int(x) for x in (args.StereoCentres).split(',')]
    if args.GenOnly:
        settings.GenOnly = True
    if args.AssumeDFTDone:
        settings.AssumeDone = True
    if args.solvent:
        settings.Solvent = args.solvent
    if args.rot5:
        settings.Rot5Cycle = True
    if args.ra is not None:
        settings.RingAtoms =\
            [int(x) for x in (args.ra).split(',')]
    if len(args.StructureFiles) == 1:
        main(args.StructureFiles[0], args.ExpNMR, 1)
    else:
        main(args.StructureFiles, args.ExpNMR, len(args.StructureFiles))
    #main()
