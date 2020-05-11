# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 13:15:47 2015

@author: ke291

Contains all of the MacroModel specific code for input generation, calculation
execution and output interpretation. Called by PyDP4.py.
"""

import os
import shutil
import sys
import subprocess
import shutil
import time
import re


def SetupMacroModel(settings):

    if settings.SCHRODINGER == '':
        SchrodEnv = os.getenv('SCHRODINGER')
        if SchrodEnv != None:
            settings.SCHRODINGER = SchrodEnv
        else:
            if os.path.exists('/usr/local/shared/schrodinger/current'):
                settings.SCHRODINGER = '/usr/local/shared/schrodinger/current'
            else:
                print('Could not find Schrodinger folder, please fill in PyDP4.py Settings with the path.')

    MacroModelInputs = []

    for f in settings.InputFiles:
        
        if settings.Rot5Cycle is True:
            if not os.path.exists(f+'rot.sdf'):
                import FiveConf
                #Generate the flipped fivemembered ring,
                #result is in '*rot.sdf' file
                FiveConf.main(f + '.sdf', settings)

        scriptdir = getScriptPath()
        cwd = os.getcwd()

        #Convert input structure to mae file
        if os.name == 'nt':
            convinp = '"' + settings.SCHRODINGER + '/utilities/sdconvert" -isd '
        else:
            convinp = settings.SCHRODINGER + '/utilities/sdconvert -isd '

        if (f[-4:] == '.sdf'):
            if not os.path.exists(f[:-4] + '.mae'):
                outp = subprocess.check_output(convinp + f + ' -omae ' + f[:-4] +
                                               '.mae', shell=True)
            MacroModelInputs.append(f[:-4] + '.mae')
        else:
            if not os.path.exists(f + '.mae'):
                outp = subprocess.check_output(convinp + f + '.sdf -omae ' + f +
                                               '.mae', shell=True)
            MacroModelInputs.append(f + '.mae')

        #Copy default com file to directory
        shutil.copyfile(settings.ScriptDir + '/default.com', cwd + '/' + f + '.com')
        #Change input and output file names in com file
        comf = open(f + '.com', 'r+')
        com = comf.readlines()
        com[0] = f + '.mae\n'
        com[1] = f + '-out.mae\n'

        #Change the molecular mechanics step count in the com file
        cycles = (str(settings.MMstepcount)).rjust(6)
        temp = list(com[7])
        temp[7:13] = list(cycles)
        com[7] = "".join(temp)
        comf.truncate(0)
        comf.seek(0)
        comf.writelines(com)
                
        #Change the forcefield in the com file
        if (settings.ForceField).lower() == "opls":
            temp = list(com[3])
            
            temp[11:13] = list('14')
            com[3] = "".join(temp)
            comf.truncate(0)
            
            comf.seek(0)
            comf.writelines(com)
        
        comf.close()
        
        if settings.Rot5Cycle is True:
            #Convert input structure to mae file
            if os.name == 'nt':
                convinp = '"' + settings.SCHRODINGER + '/utilities/sdconvert" -isd '
            else:
                convinp = settings.SCHRODINGER + '/utilities/sdconvert -isd '
            outp = subprocess.check_output(convinp + f + 'rot.sdf -omae ' + f +
                                           'rot.mae', shell=True)
            MacroModelInputs.append(f + 'rot.mae')

            #Copy default com file to directory
            shutil.copyfile(settings.ScriptDir + '/default.com', cwd + '/' + f +
                            'rot.com')
            #Change input and output file names in com file
            comf = open(f + 'rot.com', 'r+')
            com = comf.readlines()
            com[0] = f + 'rot.mae\n'
            com[1] = f + 'rot-out.mae\n'

            #Change the molecular mechanics step count in the com file
            cycles = (str(settings.MMstepcount)).rjust(6)
            temp = list(com[7])
            temp[7:13] = list(cycles)
            com[7] = "".join(temp)
            comf.truncate(0)
            comf.seek(0)
            comf.writelines(com)
            comf.close()
        print("Macromodel input for " + f + " prepared.")

    return MacroModelInputs


def RunMacroModel(MacroModelInputs, settings):
    #not args, but MacroModelInputs, numDS removed
    #Run Macromodel conformation search for all diastereomeric inputs

    MacroModelBaseNames = [x[:-4] for x in MacroModelInputs]
    MacroModelOutputs = []
    NCompleted = 0

    if shutil.which(os.path.join(settings.SCHRODINGER,'bmin')) is None & settings.SCHRODINGER != "not-installed":
        print('MacroModel.py, RunMacroModel:\n  Could not find MacroModel executable at ' +
              os.path.join(settings.SCHRODINGER,'bmin'))
        quit()

    if os.name == 'nt':
        MMPrefix = '"' + settings.SCHRODINGER + '\\bmin" '
    else:
        MMPrefix = settings.SCHRODINGER + "/bmin "

    for isomer in MacroModelBaseNames:
        if not os.path.exists(isomer + '.log'):
            print(MMPrefix + isomer)
            outp = subprocess.check_output(MMPrefix + isomer, shell=True)
        else:
            if IsMMCompleted(isomer + '.log'):
                print("Valid " + isomer + ".log exists, skipping...")
                NCompleted = NCompleted + 1
                MacroModelOutputs.append(isomer + '.log')
            else:
                print("Incomplete " + isomer + ".log exists, consider deleting it. Skipping...")
            continue

        time.sleep(60)
        while(not IsMMCompleted(isomer + '.log')):
            time.sleep(30)
        NCompleted = NCompleted + 1
        MacroModelOutputs.append(isomer + '.log')

        print("Macromodel job " + str(NCompleted) + " of " + str(len(MacroModelBaseNames)) + " completed.")

    return MacroModelOutputs


def ReadConformers(MacroModelOutputs, Isomers, settings):
    MatchingOutput = ''
    for iso in Isomers:
        for outp in MacroModelOutputs:
            if (outp[:-4] == iso.BaseName) and IsMMCompleted(outp):
                print(outp + ' is matching conformational search output for ' + iso.BaseName)
                atoms, conformers, charge, AbsEs = ReadMacromodel(iso.BaseName, settings)
                iso.Atoms = atoms
                iso.Conformers = conformers
                iso.MMCharge = charge
                iso.MMEnergies = AbsEs
    return Isomers


def ReadMacromodel(MMoutp, settings):

    conformers = []
    conformer = -1
    AbsEs = []
    ConfAbsEs = []

    atoms = []
    charge = 0
    MaeInps = []

    MaeFile = open(MMoutp + '-out.mae', 'r')
    MaeInp = MaeFile.readlines()
    MaeFile.close()

    AbsEOffsets = []
    #find conformer description blocks
    blocks = []
    DataIndexes = []
    for i in range(len(MaeInp)):
        if 'f_m_ct' in MaeInp[i]:
            blocks.append(i)
        if 'p_m_ct' in MaeInp[i]:
            blocks.append(i)

    #find absolute energy offsets
    for block in blocks:
        for i in range(block, len(MaeInp)):
            if 'mmod_Potential_Energy' in MaeInp[i]:
                AbsEOffsets.append(i-block)
                break

    #Get absolute energies for conformers
    for i in range(len(blocks)):
        for line in range(blocks[i], len(MaeInp)):
            if ':::' in MaeInp[line]:
                AbsEs.append(float(MaeInp[line+AbsEOffsets[i]]))
                break

    #Pick only the conformers in the energy window
    MinE = min(AbsEs)

    #find geometry descriptions for each block
    for i in range(len(blocks)):
        for line in (MaeInp[blocks[i]:]):
            if 'm_atom' in line:
                blocks[i] = blocks[i] + MaeInp[blocks[i]:].index(line)
                break

    #find start of atom coordinates for each block
    for i in range(len(blocks)):
        if (AbsEs[i] < MinE+settings.MaxCutoffEnergy):
            #Save the locations of atom number, xyz and charge
            DataIndex = [0, 0, 0, 0, 0]
            for offset, line in enumerate(MaeInp[blocks[i]:]):
                if 'i_m_mmod_type' in line:
                    DataIndex[0] = offset-1
                if 'r_m_x_coord' in line:
                    DataIndex[1] = offset-1
                if 'r_m_y_coord' in line:
                    DataIndex[2] = offset-1
                if 'r_m_z_coord' in line:
                    DataIndex[3] = offset-1
                if 'r_m_charge1' in line:
                    DataIndex[4] = offset-1
                if ':::' in line:
                    blocks[i] = blocks[i] + MaeInp[blocks[i]:].index(line)
                    break
            DataIndexes.append(DataIndex)
        else:
            break

    #Read the atom numbers and coordinates
    for i, block in enumerate(blocks):
        if (AbsEs[i] < MinE+settings.MaxCutoffEnergy):
            conformers.append([])
            ConfAbsEs.append(AbsEs[i])
            conformer = conformer + 1
            index = block+1
            atom = 0
            while not ':::' in MaeInp[index]:
                #Replace quoted fields with x
                line = (re.sub(r"\".*?\"", "x", MaeInp[index],
                                    flags=re.DOTALL)).split(' ')
                line = [word for word in line[:-1] if word != '']
                conformers[conformer].append([])
                if conformer == 0:
                    atoms.append(GetMacromodelSymbol(int(line[DataIndexes[i][0]])))
                    conformers[0][atom].append(line[DataIndexes[i][1]])  # add X
                    conformers[0][atom].append(line[DataIndexes[i][2]])  # add Y
                    conformers[0][atom].append(line[DataIndexes[i][3]])  # add Z
                    charge = charge + float(line[DataIndexes[i][4]])

                else:
                    conformers[conformer][atom].append(line[DataIndexes[i][1]])  # add X
                    conformers[conformer][atom].append(line[DataIndexes[i][2]])  # add Y
                    conformers[conformer][atom].append(line[DataIndexes[i][3]])  # add Z

                index = index + 1   # Move to next line
                atom = atom + 1     # Move to next atom
        else:
            break

    return atoms, conformers, int(charge), ConfAbsEs


def GetMacromodelSymbol(atomType):

    Lookup = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
              'C', 'C', 'C', 'O', 'O', 'O', ' ', 'O', ' ', 'O',
              'O', ' ', 'O', 'N', 'N', 'N', ' ', ' ', ' ', ' ',
              'N', 'N', ' ', ' ', ' ', ' ', ' ', 'N', 'N', 'N',
              'H', 'H', 'H', 'H', 'H', ' ', ' ', 'H', 'S', ' ',
              'S', 'S', 'P', 'B', 'B', 'F', 'Cl', 'Br', 'I', 'Si',
              ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
              ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
              ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
              ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'S',
              'S', 'Cl', 'B', 'F', ' ', ' ', ' ', ' ', 'S', 'S',
              ' ', ' ', 'S', 'S']

    if Lookup[atomType-1] == ' ':
        print('Unknown atom type')

    return Lookup[atomType-1]


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def IsMMCompleted(f):
    Gfile = open(f, 'r')
    outp = Gfile.readlines()
    Gfile.close()

    if os.name == 'nt':
        i = -2
    else:
        i = -3

    if "normal termination" in outp[i]:
        return True
    else:
        return False
