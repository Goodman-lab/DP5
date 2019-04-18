#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 15:56:54 2014
Rewritten during April 2019

@author: ke291

Contains all of the Gaussian specific code for input generation and calculation
execution. Called by PyDP4.py.
"""

import subprocess
import os
import time
import glob


def SetupNMRCalcs(Isomers, settings):

    jobdir = os.getcwd()

    if not os.path.exists('nmr'):
        os.mkdir('nmr')
    os.chdir('nmr')

    for iso in Isomers:
        if iso.ExtCharge > -10:
            charge = iso.ExtCharge
        else:
            charge = iso.MMCharge

        for num in range(0, len(iso.Conformers)):
            filename = iso.BaseName + 'ginp' + str(num + 1).zfill(3)

            if os.path.exists(filename + '.out'):
                if IsGausCompleted(filename + '.out'):
                    iso.NMROutputFiles.append(filename + '.out')
                    continue
                else:
                    os.remove(filename + '.out')

            WriteGausFile(filename, iso.Conformers[num], iso.Atoms, charge, settings, 'nmr')
            iso.NMRInputFiles.append(filename + '.com')

    os.chdir(jobdir)

    return Isomers


def SetupECalcs(Isomers, settings):

    jobdir = os.getcwd()

    if not os.path.exists('e'):
        os.mkdir('e')
    os.chdir('e')

    for iso in Isomers:
        if iso.ExtCharge > -10:
            charge = iso.ExtCharge
        else:
            charge = iso.MMCharge

        for num in range(0, len(iso.Conformers)):
            filename = iso.BaseName + 'ginp' + str(num + 1).zfill(3)

            if os.path.exists(filename + '.out'):
                if IsGausCompleted(filename + '.out'):
                    iso.EOutputFiles.append(filename + '.out')
                    continue
                else:
                    os.remove(filename + '.out')

            WriteGausFile(filename, iso.Conformers[num], iso.Atoms, charge, settings, 'e')
            iso.EInputFiles.append(filename + '.com')

    os.chdir(jobdir)

    return Isomers


def RunNMRCalcs(Isomers, settings):

    jobdir = os.getcwd()
    os.chdir('nmr')

    GausJobs = []

    for iso in Isomers:
        GausJobs.extend([x for x in iso.NMRInputFiles if (x[:-4] + '.out') not in iso.NMROutputFiles])

    NCompleted = 0
    Completed = []
    gausdir = os.environ['GAUSS_EXEDIR']
    GausPrefix = gausdir + "/g09 < "

    for f in GausJobs:
        time.sleep(3)
        print(GausPrefix + f + ' > ' + f[:-3] + 'out')
        outp = subprocess.check_output(GausPrefix + f + ' > ' + f[:-3] + 'out', shell=True)

        NCompleted += 1
        if IsGausCompleted(f[:-4] + '.out'):
            Completed.append(f[:-4] + '.out')
            print("Gaussian job " + str(NCompleted) + " of " + str(len(GausJobs)) + \
                " completed.")
        else:
            print("Gaussian job terminated with an error. Continuing.")

    for iso in Isomers:
        iso.NMROutputFiles.extend([x[:-4] + '.out' for x in iso.NMRInputFiles if (x[:-4] + '.out') in Completed])

    os.chdir(jobdir)

    return Isomers


def RunECalcs(Isomers, settings):

    jobdir = os.getcwd()
    os.chdir('e')

    GausJobs = []

    for iso in Isomers:
        GausJobs.extend([x for x in iso.EInputFiles if (x[:-4] + '.out') not in iso.EOutputFiles])

    NCompleted = 0
    Completed = []
    gausdir = os.environ['GAUSS_EXEDIR']
    GausPrefix = gausdir + "/g09 < "

    for f in GausJobs:
        time.sleep(3)
        print(GausPrefix + f + ' > ' + f[:-3] + 'out')
        outp = subprocess.check_output(GausPrefix + f + ' > ' + f[:-3] + 'out', shell=True)

        NCompleted += 1
        if IsGausCompleted(f[:-4] + '.out'):
            Completed.append(f[:-4] + '.out')
            print("Gaussian job " + str(NCompleted) + " of " + str(len(GausJobs)) + \
                " completed.")
        else:
            print("Gaussian job terminated with an error. Continuing.")

    for iso in Isomers:
        iso.EOutputFiles.extend([x[:-4] + '.out' for x in iso.EInputFiles if (x[:-4] + '.out') in Completed])

    os.chdir(jobdir)

    return Isomers


def WriteGausFile(Gausinp, conformer, atoms, charge, settings, type):

    f = open(Gausinp + '.com', 'w')
    if(settings.nProc > 1):
        f.write('%nprocshared=' + str(settings.nProc) + '\n')
    if settings.DFT == 'g':
        f.write('%mem=2000MB\n%chk='+Gausinp + '.chk\n')
    else:
        f.write('%mem=6000MB\n%chk='+Gausinp + '.chk\n')

    if type == 'nmr':
        f.write(NMRroute(settings))
    elif type == 'e':
        f.write(Eroute(settings))
    elif type == 'opt':
        f.write(Optroute(settings))

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


def Optroute(settings):

    route = '# ' + settings.oFunctional + '/' + settings.oBasisSet

    if (settings.oFunctional).lower() == 'm062x':
        route += ' int=ultrafine'

    route += ' Opt=(maxcycles=' + str(settings.MaxDFTOptCycles)
    if settings.CalcFC == True:
        route += ',CalcFC'
    route += ')'

    if settings.Solvent != '':
        route += ' scrf=(solvent=' + settings.Solvent + ')'

    route += '\n'


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
        #WriteGausFileOpt(f[:-8], coords,atoms,charge,settings)
        print(f[:-8] + '* deleted and new .com files written')
    if not os.path.exists('Reoptimized.log'):
        f = file('Reoptimized.log', 'w')
        f.write('\n'.join([x[:-8] for x in GoutpFiles]))
        f.close()


#Read energy from e, if not present, then o, if not present, then nmr
def ReadDFTEnergies(Isomers, settings):
    jobdir = os.getcwd()

    if 'e' in settings.Workflow:
        os.chdir('e')
    elif 'o' in settings.Workflow:
        os.chdir('opt')
    else:
        os.chdir('nmr')

    for iso in Isomers:

        if 'e' in settings.Workflow:
            GOutpFiles = iso.EOutputFiles
        elif 'o' in settings.Workflow:
            GOutpFiles = iso.OptOutputFiles
        else:
            GOutpFiles = iso.NMROutputFiles

        DFTenergies = []
        for GOutpFile in GOutpFiles:
            gausfile = open(GOutpFile, 'r')
            GOutp = gausfile.readlines()
            gausfile.close()

            for line in GOutp:
                if 'SCF Done:' in line:
                    start = line.index(') =')
                    end = line.index('A.U.')
                    energy = float(line[start + 4:end])

            iso.DFTEnergies.append(energy)

    os.chdir(jobdir)
    return Isomers


def ReadShieldings(Isomers):

    jobdir = os.getcwd()
    os.chdir('nmr')

    for iso in Isomers:

        for GOutpFile in iso.NMROutputFiles:
            gausfile = open(GOutpFile, 'r')
            GOutp = gausfile.readlines()
            gausfile.close()

            index = 0
            shieldings = []
            labels = []

            # Find the NMR shielding calculation section
            while not 'Magnetic shielding' in GOutp[index]:
                index = index + 1

            # Read shielding constants and labels
            for line in GOutp[index:]:
                if 'Isotropic' in line:
                    data = [_f for _f in line.split(' ') if _f]
                    shieldings.append(float(data[4]))
                    labels.append(data[1] + data[0])

            iso.ConformerShieldings.append(shieldings)

        iso.ShieldingLabels = labels

    os.chdir(jobdir)
    return Isomers


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
