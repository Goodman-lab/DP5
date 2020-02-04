# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 12:52:21 2015

@author: ke291

Contains all of the Tinker specific code for input generation, calculation
execution and output interpretation. Called by PyDP4.py.
"""

import os
import sys
import subprocess


def SetupTinker(settings, *args):

    print(args)
    TinkerInputs = []
    for inpf in settings.InputFiles:
        if settings.Rot5Cycle is True:
            if not os.path.exists(inpf+'rot.sdf'):
                import FiveConf
                #Generate the flipped fivemembered ring,
                #result is put in '*rot.sdf' file
                FiveConf.main(inpf + '.sdf', settings)

        #Makes sure that the name in the title line matches filename
        #sdf2tinkerxyz uses this as the output file name
        f = open(inpf + '.sdf', 'r+')
        sdf = f.readlines()
        sdf[0] = inpf + '\n'
        f.seek(0)
        f.writelines(sdf)
        f.close()

        if settings.Rot5Cycle is True:
            f = open(inpf + 'rot.sdf', 'r+')
            sdf = f.readlines()
            sdf[0] = inpf + 'rot\n'
            if inpf in sdf[-3]:
                sdf[-3] = inpf + 'rot.1\n'
            f.seek(0)
            f.writelines(sdf)
            f.close()

        scriptdir = getScriptPath()
        print(scriptdir)
        convinp = scriptdir + '/sdf2tinkerxyz -k ' + scriptdir + '/default.key <'

        outp = subprocess.check_output(convinp + inpf + '.sdf', shell=True)
        TinkerInputs.append(inpf + '.xyz')
        print("Tinker input for " + inpf + " prepared.")

        if settings.Rot5Cycle is True:
            outp = subprocess.check_output(convinp + inpf + 'rot.sdf',
                                           shell=True)
            print("Tinker input for " + inpf + "rot prepared.")

    return TinkerInputs


def RunTinker(numDS, settings, *args):
    #Run Tinker scan for all diastereomeric inputs

    NCompleted = 0

    for ds in args:
        print(settings.TinkerPath + ds + ' 0 10 20 0.00001 | tee ./' + ds + \
            '.tout')
        outp = subprocess.check_output(settings.TinkerPath + ds +
            ' 0 10 20 0.00001 | tee ./' + ds + '.tout', shell=True)
        NCompleted = NCompleted + 1
        print("Tinker job " + str(NCompleted) + " of " + str(numDS) + \
            " completed.")

        if settings.Rot5Cycle is True:
            print(settings.TinkerPath + ds + 'rot 0 10 20 0.00001 | tee ./' + \
                ds + 'rot.tout')
            outp = subprocess.check_output(settings.TinkerPath + ds +
                'rot 0 10 20 0.00001 | tee ./' + ds + 'rot.tout', shell=True)
            NCompleted = NCompleted + 1
            print("Tinker job " + str(NCompleted) + " of " + str(numDS*2) + \
                " completed.")


#Reads the relevant tinker geometry files
#v0.2 - reads seperate rot file as well
def ReadTinker(TinkerOutput, settings):

    #Get conformer energies
    ETable, charge = GetEnergiesCharge(TinkerOutput, settings)

    if settings.Rot5Cycle is True:
        #Get conformer energies for the flipped 5-membered ring
        ETableRot, charge = GetEnergiesCharge(TinkerOutput + 'rot', settings)

    #Determine which conformers we want
    MinE = 10000
    MinE = min([float(x[1]) for x in ETable])

    if settings.Rot5Cycle is True:
        MinERot = 10000
        MinERot = min([float(x[1]) for x in ETableRot])
        if MinE > MinERot:
            MinE = MinERot

    FileNums = []
    RotFileNums = []

    AcceptedEs = []

    for conf in ETable:
        if float(conf[1]) < MinE + settings.MaxCutoffEnergy:
             #Dealing with special case when nconf>100 000
            if 'Minimum' in conf[0]:
                data = conf[0].strip()
                FileNums.append(data[7:].strip())
            else:
                FileNums.append(conf[0].strip())
            AcceptedEs.append(float(conf[1]))

    if settings.Rot5Cycle is True:
        for conf in ETableRot:
            if float(conf[1]) < MinE + settings.MaxCutoffEnergy:
                RotFileNums.append(conf[0].strip())
                AcceptedEs.append(float(conf[1]))

    print("Number of accepted conformers by energies: " + str(len(AcceptedEs)))

    Files = []
    #Generate conformer filenames
    for num in FileNums:
        Files.append(TinkerOutput + '.' + num.zfill(3))
    if settings.Rot5Cycle is True:
        for num in RotFileNums:
            Files.append(TinkerOutput + 'rot.' + num.zfill(3))

    conformers = []
    conformer = 0
    atoms = []

    for f in Files:
        conformers.append([])

        atom = 0
        infile = open(f, 'r')
        inp = infile.readlines()

        for line in inp[1:]:
            data = line.split(' ')
            data = [_f for _f in data if _f]
            if conformer == 0:
                atoms.append(GetTinkerSymbol(int(data[5])))  # Add atom symbol
            conformers[conformer].append([])                # Add new atom
            conformers[conformer][atom].append(data[0])     # add atom number
            conformers[conformer][atom].append(data[2])     # add X
            conformers[conformer][atom].append(data[3])     # add Y
            conformers[conformer][atom].append(data[4])     # add Z
            atom = atom + 1     # Move to the next atom

        infile.close()
        conformer = conformer + 1   # Move to the next conformer

    return atoms, conformers, charge


# Get energies of conformers from tinker output file
def GetEnergiesCharge(TinkerOutput, settings):

    infile = open(TinkerOutput + '.tout', 'r')

    inp = infile.readlines()
    if len(inp) < 56:
        print("Tinker output " + TinkerOutput + " is corrupted, aborting.")
        quit()

    #Get the conformer energies from the file
    ETable = []
    for line in inp[13:]:
        data = line.split('  ')
        data = [_f for _f in data if _f]
        if len(data) >= 3:
            if 'Map' in data[0] and 'Minimum' in data[1]:
                ETable.append(data[-2:])
                #print data

    infile.close()
    if settings.charge is None:
        if os.path.exists(TinkerOutput+'.inchi'):
            return ETable, GetInchiCharge(TinkerOutput)
        else:
            return ETable, GetSDFCharge(TinkerOutput)
    else:
        return ETable, settings.charge


#translate Tinker atom types to element symbols for NWChem file
def GetTinkerSymbol(atomType):

    Lookup = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', ' ', 'C',
              'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
              'C', 'C', 'H', 'H', 'O', 'O', 'O', 'O', 'O', 'O',
              'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
              'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'N', 'N',
              'N', 'N', 'N', 'N', 'N', 'F', 'Cl', 'Br', 'I', 'S',
              'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S', 'S',
              'Si', 'C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'H',
              'H', 'H', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'P',
              'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
              'H', 'H', 'H', 'H', 'C', 'H', 'O', 'O', 'O', 'O',
              'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O',
              'O', 'H', 'N', 'O', 'O', 'H', 'H', 'H', 'H', 'H',
              'H', 'H', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'C',
              'C', 'N', 'N', 'N', 'N', 'N', 'N', 'S', 'N', 'N',
              'N', 'N', 'N', 'O', 'H', 'O', 'H', 'N', 'N', 'N',
              'N', 'N', 'C', 'C', 'N', 'O', 'C', 'N', 'N', 'C',
              'C', 'N', 'N', 'N', 'N', 'N', 'O', 'H', 'H', 'H',
              'S', 'S', 'S', 'S', 'S', 'S', 'S', 'P', 'N', 'Cl',
              'C', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
              'Fe', 'Fe', 'F', 'Cl', 'Br', 'Li', 'Na', 'K', 'Zn', 'Zn',
              'Ca', 'Cu', 'Cu', 'Mg']

    if Lookup[atomType-1] == ' ':
        print('Unknown atom type')

    return Lookup[atomType-1]


def GetInchiCharge(inchifile):

    infile = open(inchifile + '.inchi', 'r')
    inp = infile.readlines()
    infile.close()

    ChargeFound = False
    #Get inchi layers
    layers = inp[0].split('/')
    for l in layers[1:]:
        if 'q' in l:
            charge = int(l[1:])
            ChargeFound = True
            break
        if 'p' in l:
            charge = int(l[1:])
            ChargeFound = True
            break
    if not ChargeFound:
        charge = 0

    return charge


def GetSDFCharge(sdf):
    import openbabel

    obconversion = openbabel.OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = openbabel.OBMol()
    obconversion.ReadFile(obmol, sdf+'.sdf')

    return obmol.GetTotalCharge()


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
