# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 13:15:47 2015

@author: ke291

Contains all of the MacroModel specific code for input generation, calculation
execution and output interpretation. Called by PyDP4.py.
"""

import os
import sys
import subprocess
import shutil
import time
import re


def SetupMacromodel(numDS, settings, *args):

    for f in args[0]:
        
        if settings.Rot5Cycle is True:
            if not os.path.exists(f+'rot.sdf'):
                import FiveConf
                #Generate the flipped fivemembered ring,
                #result is in '*rot.sdf' file
                FiveConf.main(f + '.sdf', settings)

        scriptdir = getScriptPath()
        cwd = os.getcwd()

        #Convert input structure to mae file
        convinp = settings.SCHRODINGER + '/utilities/sdconvert -isd '
        outp = subprocess.check_output(convinp + f + '.sdf -omae ' + f +
                                       '.mae', shell=True)

        #Copy default com file to directory
        shutil.copyfile(scriptdir + '/default.com', cwd + '/' + f + '.com')
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
            convinp = settings.SCHRODINGER + '/utilities/sdconvert -isd '
            outp = subprocess.check_output(convinp + f + 'rot.sdf -omae ' + f
                                           + 'rot.mae', shell=True)

            #Copy default com file to directory
            shutil.copyfile(scriptdir + '/default.com', cwd + '/' + f +
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
        print "Macromodel input for " + f + " prepared."


def RunMacromodel(numDS, settings, *args):
    #Run Macromodel conformation search for all diastereomeric inputs
    NCompleted = 0

    for ds in args[0]:
        if not os.path.exists(ds+'.log'):
            print settings.SCHRODINGER + '/bmin ' + ds
            outp = subprocess.check_output(settings.SCHRODINGER + '/bmin ' +
                                           ds, shell=True)
        else:
            print ds + ".log exists, skipping"
            continue

        time.sleep(60)
        while(not IsMMCompleted(ds + '.log')):
            time.sleep(30)
        NCompleted = NCompleted + 1

        if settings.Rot5Cycle is True:
            print "Macromodel job " + str(NCompleted) + " of " + str(numDS*2)\
                + " completed."

            print settings.SCHRODINGER + '/bmin ' + ds + 'rot'
            outp = subprocess.check_output(settings.SCHRODINGER + '/bmin '
                                           + ds+'rot', shell=True)
            time.sleep(60)
            while(not IsMMCompleted(ds + 'rot.log')):
                time.sleep(30)
            NCompleted = NCompleted + 1
            print "Macromodel job " + str(NCompleted) + " of " + str(numDS*2)\
                + " completed."
        else:
            print "Macromodel job " + str(NCompleted) + " of " + str(numDS) +\
                " completed."


def ReadMacromodel(MMoutp, settings):

    conformers = []
    conformer = -1
    AbsEs = []

    atoms = []
    charge = 0
    MaeInps = []

    MaeFile = file(MMoutp + '-out.mae', 'r')
    MaeInps.append(MaeFile.readlines())
    MaeFile.close()

    if settings.Rot5Cycle is True:
        MaeFile = file(MMoutp + 'rot-out.mae', 'r')
        MaeInps.append(MaeFile.readlines())
        MaeFile.close()

    for MaeInp in MaeInps:
        index = 0
        AbsEOffsets = []
        #find conformer description blocks
        blocks = []
        for i in range(0, len(MaeInp)):
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
        for i in range(0, len(blocks)):
            for line in range(blocks[i], len(MaeInp)):
                if ':::' in MaeInp[line]:
                    AbsEs.append(float(MaeInp[line+AbsEOffsets[i]]))
                    break

        #find geometry descriptions for each block
        for i in range(0, len(blocks)):
            for line in (MaeInp[blocks[i]:]):
                if 'm_atom' in line:
                    blocks[i] = blocks[i] + MaeInp[blocks[i]:].index(line)
                    break

        #find start of atom coordinates for each block
        for i in range(0, len(blocks)):
            for line in (MaeInp[blocks[i]:]):
                if ':::' in line:
                    blocks[i] = blocks[i] + MaeInp[blocks[i]:].index(line)
                    break

        #Read the atom numbers and coordinates
        for block in blocks:
            conformers.append([])
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
                    atoms.append(GetMacromodelSymbol(int(line[1])))
                    conformers[0][atom].append(line[0])  # add atom number
                    conformers[0][atom].append(line[2])  # add X
                    conformers[0][atom].append(line[3])  # add Y
                    conformers[0][atom].append(line[4])  # add Z
                    charge = charge + int(line[16])
                    
                else:
                    if blocks.index(block) == 0:
                        conformers[conformer][atom].append(line[0])  # add atom number
                        conformers[conformer][atom].append(line[2])  # add X
                        conformers[conformer][atom].append(line[3])  # add Y
                        conformers[conformer][atom].append(line[4])  # add Z
                    else:
                        conformers[conformer][atom].append(line[0])  # add atom number
                        conformers[conformer][atom].append(line[1])  # add X
                        conformers[conformer][atom].append(line[2])  # add Y
                        conformers[conformer][atom].append(line[3])  # add Z

                index = index + 1   # Move to next line
                atom = atom + 1     # Move to next atom
    #Pick only the conformers in the energy window
    MinE = min(AbsEs)

    conformers2 = []
    for i in range(0, len(conformers)):
        if AbsEs[i] < MinE+settings.MaxCutoffEnergy:
            conformers2.append(conformers[i])

    return atoms, conformers2, charge


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
        print 'Unknown atom type'

    return Lookup[atomType-1]


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def IsMMCompleted(f):
    Gfile = open(f, 'r')
    outp = Gfile.readlines()
    Gfile.close()

    if "normal termination" in outp[-3]:
        return True
    else:
        return False
