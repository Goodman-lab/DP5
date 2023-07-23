# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 12:52:21 2015

@author: ke291

Contains all of the Tinker specific code for input generation, calculation
execution and output interpretation. Called by PyDP4.py.
"""

import os
import shlex
import shutil
import sys
import subprocess
import PyDP4

# Please modify the line below to give the path to the TINKER v8.x top level folder
# This folder should contain bin/scan and params/mmff.prm for the process to work

def SetupTinker(settings):

    TinkerInputs = []
    for inpfi in settings.InputFiles:

        inpf = inpfi.split('.')[0]

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
        #convinp = scriptdir + '/sdf2tinkerxyz -k ' + scriptdir + '/default.key <'
        #outp = subprocess.check_output(convinp + inpf + '.sdf', shell=True)

        import sdftinkerxyzpy

        convinp = sdftinkerxyzpy.main(inpf)

        #f = open(inpf + '.key', 'w+')
        #key = f.readlines()
        #key[2] = 'PARAMETERS        ' + settings.TinkerPath + 'params/mmff.prm\n'
        #f.seek(0)
        #f.writelines(key)
        #f.close()

        TinkerInputs.append(inpf)
        print("Tinker input for " + inpf + " prepared.")

        if settings.Rot5Cycle is True:
            outp = subprocess.check_output(convinp + inpf + 'rot.sdf',
                                           shell=True)
            print("Tinker input for " + inpf + "rot prepared.")

    return TinkerInputs

def RunTinker(TinkerInputs, settings):
    #Run Tinker scan for all diastereomeric inputs
    TinkerOutputs = []
    TinkerPrefix = os.path.join(settings.TinkerPath, 'bin', 'scan')
    if shutil.which(TinkerPrefix) is None:
        print('Tinker.py, RunTinker:\n  Could not find Tinker scan executable at ' + TinkerPrefix)
        quit()

    NCompleted = 0

    for isomer in TinkerInputs:
        if os.path.exists(isomer + '.tout') and os.path.exists(isomer + '.arc'):
            print('Output files for ' + isomer + ' already exist, skipping.')
            TinkerOutputs.append(isomer)
            continue

        print(settings.TinkerPath + '/bin/scan ' + isomer +' '+settings.TinkerPath + '/params/mmff.prm 0 10 20 0.00001 | tee ./' + isomer + \
            '.tout')
        outp = subprocess.check_output(settings.TinkerPath + '/bin/scan ' + isomer +' '+settings.TinkerPath+'/params/mmff.prm'+
            ' 0 10 20 0.00001 | tee ./' + isomer + '.tout', shell=True)
        NCompleted = NCompleted + 1
        TinkerOutputs.append(isomer)
        print("Tinker job " + str(NCompleted) + " of " + str(len(TinkerInputs)) + \
            " completed.")

        if settings.Rot5Cycle is True:
            print(settings.TinkerPath + '/bin/scan ' + isomer + 'rot '+settings.TinkerPath+'/params/mmff.prm 0 10 20 0.00001 | tee ./' + \
                isomer + 'rot.tout')
            outp = subprocess.check_output(settings.TinkerPath + '/bin/scan ' + isomer +
                'rot '+settings.TinkerPath+'/params/mmff.prm 0 10 20 0.00001 | tee ./' + isomer + 'rot.tout', shell=True)
            NCompleted = NCompleted + 1

    return TinkerOutputs


def ReadConformers(TinkerOutputs, Isomers, settings):
    atypes, anums = ExtractAtomTypes(settings)

    for iso in Isomers:
        for outp in TinkerOutputs:
            if (outp == iso.BaseName):
                print(outp + ' is matching conformational search output for ' + iso.BaseName)
                atoms, conformers, charge, AbsEs = ReadTinker(outp, settings, atypes, anums)
                iso.Atoms = atoms
                iso.Conformers = conformers
                iso.MMCharge = charge
                iso.MMEnergies = AbsEs
            else:
                print(outp + ' ' + iso.BaseName + ' Nope')
    return Isomers


#Reads force field parameter file to understand atom notation in the output
def ExtractAtomTypes(settings):
    # open settings.TinkerPath + 'params/mmff.prm
    atomtypes = []
    atomnums = []
    with open(settings.TinkerPath + '/params/mmff.prm', 'r') as f:
        for line in f:
            if line.split(' ')[0] == 'atom':
                data = shlex.split(line, posix=False)
                atomtypes.append(data[3])
                atomnums.append(int(data[-3]))
    return atomtypes, atomnums

#Reads the relevant tinker geometry files
#v0.2 - reads seperate rot file as well
def ReadTinker(TinkerOutput, settings, atypes, anums):

    print('Reading ' + TinkerOutput)
    #Get conformer energies
    energies, charge = GetEnergiesCharge(TinkerOutput, settings)

    if settings.Rot5Cycle is True:
        #Get conformer energies for the flipped 5-membered ring
        energiesrot, charge = GetEnergiesCharge(TinkerOutput + 'rot', settings)
        energies = energies + energiesrot
    #Determine which conformers we want
    MinE = 10000
    MinE = min(energies)

    if settings.Rot5Cycle is True:
        MinERot = 10000
        MinERot = min(energiesrot)
        if MinE > MinERot:
            MinE = MinERot

    atoms, conformers = ReadArc(TinkerOutput, atypes, anums)
    if settings.Rot5Cycle is True:
        # Get conformer energies for the flipped 5-membered ring
        atoms, conformersrot = ReadArc(TinkerOutput + 'rot', atypes, anums)
        conformers = conformers + conformersrot

    filtered = []
    AbsEs = []

    for energy, conformer in zip(energies, conformers):
        if energy < MinE + settings.MaxCutoffEnergy:
            AbsEs.append(energy)
            filtered.append(conformer)

    print("Number of accepted conformers by energies: " + str(len(filtered)))

    return atoms, filtered, charge, AbsEs


# Reads all conformers from the arc file
def ReadArc(f, atypes, anums):
    conffile = open(f + '.arc', 'r')
    confdata = conffile.readlines()
    conffile.close()
    #output data: conformers - list of x,y,z lists, atoms - list of atoms
    conformers = []
    atoms = []
    atypes = [x[:3] for x in atypes]

    for line in confdata:
        data = [_f for _f in line.split('  ') if _f]
        if len(data) < 3:
            conformers.append([])
        else:
            if len(conformers) == 1:
                anum = anums[atypes.index(data[1][:3])]
                atoms.append(GetAtomSymbol(anum))
            conformers[-1].append([x for x in data[2:5]])

    return atoms, conformers


# Get energies of conformers from tinker output file
def GetEnergiesCharge(TinkerOutput, settings):

    infile = open(TinkerOutput + '.tout', 'r')

    inp = infile.readlines()
    if len(inp) < 13:
        print("Tinker output " + TinkerOutput + " is corrupted, aborting.")
        quit()

    #Get the conformer energies from the file
    energies = []
    for line in inp[13:]:
        data = line[:-1].split('  ')
        data = [_f for _f in data if _f]
        if len(data) >= 3:
            if 'Map' in data[0] and 'Minimum' in data[1]:
                energies.append(float(data[-1]))
                #print data

    infile.close()
    if settings.charge is None:
        if os.path.exists(TinkerOutput+'.inchi'):
            return energies, GetInchiCharge(TinkerOutput)
        else:
            return energies, GetSDFCharge(TinkerOutput)
    else:
        return energies, settings.charge


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
