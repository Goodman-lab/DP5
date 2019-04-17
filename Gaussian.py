#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 15:56:54 2014

@author: ke291

Contains all of the Gaussian specific code for input generation and calculation
execution. Called by PyDP4.py.
"""

import nmrPredictGaus

import subprocess
import os
import time
import glob


def SetupNMRCalcs(Isomers, settings):

    jobdir = os.getcwd()

    if not os.path.exists('./nmr'):
        os.mkdir('./nmr')
    os.chdir('./nmr')

    for iso in Isomers:
        if iso.ExtCharge > -10:
            charge = iso.ExtCharge
        else:
            charge = iso.MMCharge

        for num in range(0, len(iso.Conformers)):
            filename = iso.BaseName + 'ginp' + str(num + 1).zfill(3)

            if os.path.exists(filename + '.out') and IsGausCompleted(filename + '.out'):
                continue

            WriteGausFile(filename, iso.Conformers[num], iso.Atoms, charge, settings)

    os.chdir(jobdir)


def WriteGausFile(Gausinp, conformer, atoms, charge, settings):

    f = open(Gausinp + '.com', 'w')
    if(settings.nProc > 1):
        f.write('%nprocshared=' + str(settings.nProc) + '\n')
    if settings.DFT == 'g':
        f.write('%mem=2000MB\n%chk='+Gausinp + '.chk\n')
    else:
        f.write('%mem=6000MB\n%chk='+Gausinp + '.chk\n')

    f.write(NMRroute(settings))
    f.write('\n'+Gausinp+'\n\n')
    f.write(str(charge) + ' 1\n')

    natom = 0

    for atom in conformer:
        f.write(atoms[natom] + '  ' + atom[0] + '  ' + atom[1] + '  ' +
                atom[2] + '\n')
        natom = natom + 1
    f.write('\n')

    f.close()


def NMRroute(settings):

    route = '# ' + settings.nFunctional + '/' + settings.nBasisSet
    if (settings.nFunctional).lower() == 'm062x':
        route += ' int=ultrafine'

    route += ' nmr=giao'

    if settings.Solvent != '':
        route += ' scrf=(solvent=' + settings.Solvent + ')'

    route += '\n'

    return route


def Eroute(settings):

    route = '# ' + settings.eFunctional + '/' + settings.eBasisSet
    if (settings.eFunctional).lower() == 'm062x':
        route += ' int=ultrafine'

    if settings.Solvent != '':
        route += ' scrf=(solvent=' + settings.Solvent + ')'

    route += '\n'

    return route


def WriteGausFileOpt(Gausinp, conformer, atoms, charge, settings):

    #write the initial DFT geometry optimisation input file first
    f1 = file(Gausinp + 'a.com', 'w')
    if(settings.nProc > 1):
        f1.write('%nprocshared=' + str(settings.nProc) + '\n')

    if settings.DFT != 'd':
        f1.write('%mem=6000MB\n%chk='+Gausinp + '.chk\n')
    else:
        fullscrfolder = settings.DarwinScrDir + settings.StartTime + settings.Title + '/'
        f1.write('%mem=6000MB\n%chk=' + fullscrfolder + Gausinp + '.chk\n')
    
    if settings.DFTOpt:
        if settings.Solvent != '':
            f1.write('# b3lyp/6-31g(d,p) Opt=(maxcycles=' + str(settings.MaxDFTOptCycles) + ') scrf=(solvent=' +
                     settings.Solvent+')\n')
        else:
            f1.write('# b3lyp/6-31g(d,p) Opt=(maxcycles=' + str(settings.MaxDFTOptCycles) + ')\n')
    elif settings.PM6Opt:
        if settings.Solvent != '':
            f1.write('# pm6 Opt=(maxcycles=' + str(settings.MaxDFTOptCycles) + ') scrf=(solvent=' +
                     settings.Solvent+')\n')
        else:
            f1.write('# pm6 Opt=(maxcycles=' + str(settings.MaxDFTOptCycles) + ')\n')
    elif settings.HFOpt:
        if settings.Solvent != '':
            f1.write('# rhf/6-31g(d,p) Opt=(maxcycles=' + str(settings.MaxDFTOptCycles) + ') scrf=(solvent=' +
                     settings.Solvent+')\n')
        else:
            f1.write('# rhf/6-31g(d,p) Opt=(maxcycles=' + str(settings.MaxDFTOptCycles) + ')\n')
    elif settings.M06Opt:
        if settings.Solvent != '':
            f1.write('# m062x/6-31g(d,p) Opt=(maxcycles=' + str(settings.MaxDFTOptCycles) + ') scrf=(solvent=' +
                     settings.Solvent+')\n')
        else:
            f1.write('# m062x/6-31g(d,p) Opt=(maxcycles=' + str(settings.MaxDFTOptCycles) + ')\n')
            

    f1.write('\n'+Gausinp+'\n\n')
    f1.write(str(charge) + ' 1\n')

    natom = 0

    for atom in conformer:
        f1.write(atoms[natom] + '  ' + atom[1] + '  ' + atom[2] + '  ' +
                 atom[3] + '\n')
        natom = natom + 1
    f1.write('\n')
    f1.close()

    #Write the nmr prediction input file,
    #using the geometry from checkpoint file
    f2 = file(Gausinp + 'b.com', 'w')
    if(settings.nProc > 1):
        f2.write('%nprocshared=' + str(settings.nProc) + '\n')

    if settings.DFT != 'd':
        f2.write('%mem=6000MB\n%chk=' + Gausinp + '.chk\n')
    else:
        fullscrfolder = settings.DarwinScrDir + settings.StartTime + settings.Title + '/'
        f2.write('%mem=6000MB\n%chk=' + fullscrfolder + Gausinp + '.chk\n')

    if (settings.Functional).lower() == 'wp04':
        func1 = '# blyp/'
        func2 = ' iop(3/76=1000001189,3/77=0961409999,3/78=0000109999)' + \
            ' geom=checkpoint nmr='
        
    elif (settings.Functional).lower() == 'm062x':
        func1 = '# m062x/'
        func2 = ' int=ultrafine geom=checkpoint nmr='
        
    else:
        func1 = '# ' + settings.Functional + '/'
        func2 = ' geom=checkpoint nmr='
            
    if (settings.BasisSet).lower()[:3] == 'pcs':
        basis1 = 'gen'
    else:
        basis1 = settings.BasisSet
            
    CompSettings = func1 + basis1 + func2
    if settings.jJ or settings.jFC:
        CompSettings += '(giao,spinspin,mixed)'
    else:
        CompSettings += 'giao'
    
    if settings.Solvent != '':
        CompSettings += ' scrf=(solvent=' + settings.Solvent+')\n'
    else:
        CompSettings += '\n'
    f2.write(CompSettings)
    f2.write('\n'+Gausinp+'\n\n')
    f2.write(str(charge) + ' 1\n')
    f2.write('\n')
    if (settings.BasisSet).lower()[:3] == 'pcs':
        basisfile = file(settings.ScriptDir + '/' + 
                         (settings.BasisSet).lower(), 'r')
        inp = basisfile.readlines()
        basisfile.close()
        for line in inp:
            f2.write(line)
        f2.write('\n')
    f2.close()


def GetFiles2Run(inpfiles, settings):
    #Get the names of all relevant input files
    GinpFiles = []
    for filename in inpfiles:
        if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
            and (not settings.M06Opt):
            GinpFiles = GinpFiles + glob.glob(filename + 'ginp???.com')
        else:
            GinpFiles = GinpFiles + glob.glob(filename + 'ginp???a.com')
    Files2Run = []

    #for every input file check that there is a completed output file,
    #delete the incomplete outputs and add the inputs to be done to Files2Run
    for filename in GinpFiles:
        if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
            and (not settings.M06Opt):
            if not os.path.exists(filename[:-3]+'out'):
                Files2Run.append(filename)
            else:
                if IsGausCompleted(filename[:-3] + 'out'):
                    #print filename[:-3]+'out already exists'
                    continue
                else:
                    os.remove(filename[:-3] + 'out')
                    Files2Run.append(filename)
        else:
            if not os.path.exists(filename[:-5]+'.out'):
                Files2Run.append(filename)
            else:
                if IsGausCompleted(filename[:-5] + '.out'):
                    #print filename[:-3]+'out already exists'
                    continue
                else:
                    os.remove(filename[:-5] + '.out')
                    Files2Run.append(filename)

    return Files2Run


def IsGausCompleted(f):
    Gfile = open(f, 'r')
    outp = Gfile.readlines()
    Gfile.close()
    if len(outp) < 10:
        return False
    if "Normal termination" in outp[-1]:
        return True
    else:
        return False


def RunLocally(GausJobs, settings):
    NCompleted = 0
    gausdir = os.environ['GAUSS_EXEDIR']
    GausPrefix = gausdir + "/g09 < "

    for f in GausJobs:
        time.sleep(3)
        print(GausPrefix + f + ' > ' + f[:-3] + 'out')
        outp = subprocess.check_output(GausPrefix + f + ' > ' + f[:-3] + 'out', shell=True)
        NCompleted += 1
        print("Gaussian job " + str(NCompleted) + " of " + str(len(GausJobs)) + \
            " completed.")


def CheckConvergence(inpfiles):
    #FilesRun - list of files of the form input.com
    #we are looking for the corresponding optimization output files
    #in the form inputtemp.out
    GoutpFiles = []
    for filename in inpfiles:
        GoutpFiles = GoutpFiles + glob.glob(filename + 'ginp???temp.out')
    Nunconverged = 0
    unconverged = []
    for outfile in GoutpFiles:
        f = open(outfile, 'r')
        ginp = '\n'.join(f.readlines())
        f.close()
        if not 'Stationary point found' in ginp:
            Nunconverged += 1
            unconverged.append(outfile)
    return len(GoutpFiles), Nunconverged, unconverged


def ResubGeOpt(GoutpFiles, settings):
    for f in GoutpFiles:
        atoms, coords, charge = ReadGeometry(f[:-8]+'.out')
        for remf in glob.glob(f[:-8] + '*'):
            os.remove(remf)
        WriteGausFileOpt(f[:-8], coords,atoms,charge,settings)
        print(f[:-8] + '* deleted and new .com files written')
    if not os.path.exists('Reoptimized.log'):
        f = file('Reoptimized.log', 'w')
        f.write('\n'.join([x[:-8] for x in GoutpFiles]))
        f.close()


def ReadGeometry(GOutpFile):

    gausfile = open(GOutpFile, 'r')
    GOutp = gausfile.readlines()

    print(GOutpFile)
    print(len(GOutp))

    index = 0
    atoms = []
    coords = []

    #Find the geometry section and charge section
    for index in range(len(GOutp)):
        if 'Charge =' in GOutp[index]:
            chindex = index
        if 'Redundant internal' in GOutp[index]:
            gindex = index + 1

    #Read shielding constants and labels
    for line in GOutp[gindex:]:
        if 'Recover connectivity' in line:
            break
        else:
            data = [_f for _f in line[:-1].split(',') if _f]
            if data[0][-2:].isalpha():
                atoms.append(data[0][-2:])
            else:
                atoms.append(data[0][-1])
            coords.append(data[1:])
            
    line = GOutp[chindex].split('Charge =  ')
    line = line[1].split(' Multiplicity = ')
    charge = int(line[0])
    gausfile.close()

    return atoms, coords, charge


def ReadTempGeometry(GOutpFile):
    gausfile = open(GOutpFile, 'r')
    GOutp = gausfile.readlines()
    gausfile.close()

    if len(GOutp) < 80:
        return [], [], 0

    index = 0
    atoms = []
    coords = []
    gindex = -1
    chindex = -1

    # Find the geometry section and charge section
    for index in range(len(GOutp)):
        if 'Charge =' in GOutp[index]:
            chindex = index
        if ('Input orientation:' in GOutp[index]) or ("Standard orientation:" in GOutp[index]):
            gindex = index + 5

    if (gindex < 0) or (gindex < 0):
        return [], [], 0

    # Read shielding constants and labels
    for line in GOutp[gindex:]:
        if '--------------' in line:
            break
        else:
            data = [_f for _f in line[:-1].split(' ') if _f]
            atoms.append(GetAtomSymbol(int(data[1])))
            coords.append(data[2:])

    line = GOutp[chindex].split('Charge =  ')
    line = line[1].split(' Multiplicity = ')
    charge = int(line[0])

    return atoms, coords, charge


def GetAtomSymbol(AtomNum):
    Lookup = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', \
              'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', \
              'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', \
              'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', \
              'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', \
              'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', \
              'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']

    if AtomNum > 0 and AtomNum < len(Lookup):
        return Lookup[AtomNum-1]
    else:
        print("No such element with atomic number " + str(AtomNum))
        return 0


PTable = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', \
          'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', \
          'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', \
          'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', \
          'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', \
          'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', \
          'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']


def GetAtomNum(AtomSymbol):

    if AtomSymbol in PTable:
        return PTable.index(AtomSymbol)+1
    else:
        print("No such element with symbol " + str(AtomSymbol))
        return 0


def RunNMRPredict(numDS, settings, *args):

    GausNames = []
    NTaut = []

    for val in range(0, numDS):
        NTaut.append(args[val*2])
        GausNames.append(args[val*2+1])

    RelEs = []
    populations = []
    BoltzmannShieldings = []
    BoltzmannJs = []
    BoltzmannFCs = []
    SigConfs = []

    print(GausNames)
    print(NTaut)
    #This loop runs nmrPredict for each diastereomer and collects
    #the outputs    
    for i, isomer in enumerate(GausNames):

        GausFiles = glob.glob(isomer + 'ginp*.out')
        GausFiles = [x[:-4] for x in GausFiles]
        GausFiles = [x for x in GausFiles if 'temp' not in x]
        
        #Runs nmrPredictGaus Name001, ... and collects output
        Es, Pops, ls, BSs, Jls, BFCs, BJs, SCs = nmrPredictGaus.main(settings,
                                                                *GausFiles)
        if i == 0:
            labels = ls
            
        RelEs.append(Es)
        populations.append(Pops)
        BoltzmannShieldings.append(BSs)
        BoltzmannFCs.append(BFCs)
        BoltzmannJs.append(BJs)
        SigConfs.append(SCs)

    return (RelEs, populations, labels, BoltzmannShieldings, Jls, BoltzmannFCs,
            BoltzmannJs, SigConfs, NTaut)


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
