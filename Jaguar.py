# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 14:26:32 2015

@author: kristaps
"""
import Tinker
import MacroModel
import nmrPredictJag

import subprocess
import glob
import time
import sys
import os
import pyximport
pyximport.install()
import ConfPrune

"""
main function that runs the Jaguar, checks for when it's done and
submits the result to the NMRDP4 script for data extraction, processing
and submission to DP4 java file
"""

def SetupJaguar(MMoutp, Jaginp, numDigits, settings, adjRMSDcutoff):

    if settings.MMTinker:
        #Reads conformer geometry, energies and atom labels from Tinker output
        (atoms, conformers, charge) = Tinker.ReadTinker(MMoutp, settings)
    else:
        (atoms, conformers, charge) = MacroModel.ReadMacromodel(MMoutp,
                                                                settings)
    if settings.charge is not None:
        charge = settings.charge

    #Prune similar conformations, if the number exceeds the limit
    if len(conformers) > settings.PerStructConfLimit and settings.ConfPrune:
        pruned = ConfPrune.RMSDPrune(conformers, atoms, adjRMSDcutoff)
        actualRMSDcutoff = adjRMSDcutoff
        if len(pruned) > settings.PerStructConfLimit and settings.StrictConfLimit:
            pruned, actualRMSDcutoff = ConfPrune.StrictRMSDPrune(conformers, atoms, adjRMSDcutoff,
                                               settings.PerStructConfLimit)
    else:
        pruned = conformers
        actualRMSDcutoff = adjRMSDcutoff

    if settings.ConfPrune:
        print str(len(conformers) - len(pruned)) +\
            " or " + "{:.1f}".format(100*(len(conformers) - len(pruned)) /
            len(conformers))+"% of conformations have been pruned based on " +\
            str(actualRMSDcutoff) + " angstrom cutoff"
    
    for num in range(0, len(pruned)):
        filename = Jaginp+str(num+1).zfill(3)
        WriteJagFile(filename, pruned[num], atoms, charge, settings)
 
    print str(len(pruned)) + " .in files written"


#Adjust the RMSD cutoff to keep the conformation numbers reasonable
def AdaptiveRMSD(MMoutp, settings):

    if settings.MMTinker:
        #Reads conformer geometry, energies and atom labels from Tinker output
        (atoms, conformers, charge) = Tinker.ReadTinker(MMoutp, settings)
    else:
        (atoms, conformers, charge) = MacroModel.ReadMacromodel(MMoutp,
                                                                settings)

    return ConfPrune.AdaptRMSDPrune(conformers, atoms,
                                    settings.InitialRMSDcutoff,
                                    settings.PerStructConfLimit)



def WriteJagFile(Jaginp, conformer, atoms, charge, settings):

    f = file(Jaginp + '.in', 'w')
    
    f.write('&gen\n')
    if settings.Solvent != '':
        f.write('solvent=' + settings.Solvent+'\nisolv=2\n')
        f.write('pbf_radii=2\n')
    
    basis = settings.BasisSet
    if basis.lower() == '6-31g(d,p)':
        basis = '6-31G**'
    elif basis.lower() == '6-311g(d)':
        basis = '6-311G*'
    
    f.write('basis=' + basis + '\nnmr=1\n')
    f.write('dftname=' + settings.Functional.upper() + '\n&\n')
            
    
    #f.write('entry_name: '+Jaginp+'\n&zmat\n')
    f.write('&zmat\n')
    #f.write(str(charge) + ' 1\n')

    #natom = 0

    for natom, atom in enumerate(conformer):
        fxyz = [float(x) for x in atom[1:]]
        t = "{: .13f}"
        f.write(atoms[natom] + str(natom+1).ljust(6) + t.format(fxyz[0])\
                 + '  ' + t.format(fxyz[1]) + '  ' + t.format(fxyz[2]) + '\n')
        #natom = natom + 1
    f.write('&')
    
    f.close()


def GetFiles2Run(inpfiles, settings):
    #Get the names of all relevant input files
    
    JinpFiles = []
    for filename in inpfiles:
        JinpFiles = JinpFiles + glob.glob(filename + 'jinp???.in')
    print JinpFiles
    Files2Run = []

    #for every input file check that there is a completed output file,
    #delete the incomplete outputs and add the inputs to be done to Files2Run
    for filename in JinpFiles:
        if not os.path.exists(filename[:-2]+'out'):
            Files2Run.append(filename)
        else:
            if IsJagCompleted(filename[:-2] + 'out'):
                #print filename[:-3]+'out already exists'
                continue
            else:
                os.remove(filename[:-2] + 'out')
                Files2Run.append(filename)

    return Files2Run


def RunJaguar(JagFiles, settings):
    NCompleted = 0
    JagPrefix = settings.SCHRODINGER + "/jaguar run -WAIT "

    for f in JagFiles:
        time.sleep(3)
        print JagPrefix + f
        outp = subprocess.check_output(JagPrefix + f, shell=True)
        NCompleted += 1
        print "Jaguar job " + str(NCompleted) + " of " + str(len(JagFiles)) + \
            " completed."


def IsJagCompleted(f):
    Jfile = open(f, 'r')
    outp = Jfile.readlines()
    Jfile.close()
    if len(outp) < 10:
        return False
    if "completed on" in outp[-1]:
        return True
    else:
        return False


def RunNMRPredict(numDS, settings, *args):

    JagNames = []
    NTaut = []

    for val in range(0, numDS):
        NTaut.append(args[val*2])
        JagNames.append(args[val*2+1])

    RelEs = []
    populations = []
    BoltzmannShieldings = []
    SigConfs = []

    print JagNames
    print NTaut
    #This loop runs nmrPredict for each diastereomer and collects
    #the outputs    
    for i, isomer in enumerate(JagNames):
        print isomer
        JagFiles = glob.glob(isomer + 'jinp*.out')
        JagFiles = [x[:-4] for x in JagFiles]

        #Runs nmrPredictGaus Name001, ... and collects output
        Es, Pops, ls, BSs, SCs = nmrPredictJag.main(settings, *JagFiles)
        if i == 0:
            labels = ls
            
        RelEs.append(Es)
        populations.append(Pops)
        BoltzmannShieldings.append(BSs)
        SigConfs.append(SCs)

    return (RelEs, populations, labels, BoltzmannShieldings, SigConfs, NTaut)