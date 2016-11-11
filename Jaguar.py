# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 14:26:32 2015

@author: kristaps
"""
import Tinker
import MacroModel

import subprocess
import glob
import time
import sys
import pyximport
pyximport.install()
import ConfPrune

"""
main function that runs the Jaguar, checks for when it's done and
submits the result to the NMRDP4 script for data extraction, processing
and submission to DP4 java file
"""

def SetupJaguar(MMoutp, Gausinp, numDigits, settings, adjRMSDcutoff):

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
    
    if not settings.PM7Opt:
        for num in range(0, len(pruned)):
            filename = Gausinp+str(num+1).zfill(3)
            if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
                and (not settings.M06Opt):
                WriteGausFile(filename, pruned[num], atoms, charge, settings)
            else:
                WriteGausFileOpt(filename, pruned[num], atoms, charge, settings)
    else:
        for num in range(0, len(pruned)):
            filename = Gausinp+str(num+1).zfill(3)
            conformer = PM7opt(filename, pruned[num], atoms, charge, settings)
            print "Conformer " + str(num+1) + " of " + str(len(pruned)) + \
            " has been preoptimized at PM7 level"
            WriteGausFile(filename, conformer, atoms, charge, settings)

    print str(len(pruned)) + " .com files written"


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



def WriteJaguarFile(Gausinp, conformer, atoms, charge, settings):

    f = file(Gausinp + '.com', 'w')
    if(settings.nProc > 1):
        f.write('%nprocshared=' + str(settings.nProc) + '\n')
    f.write('%mem=6000MB\n%chk='+Gausinp + '.chk\n')
    
    if (settings.Functional).lower() == 'wp04':
        func1 = '# blyp/'
        func2 = ' iop(3/76=1000001189,3/77=0961409999,3/78=0000109999)' + \
            ' nmr='
        
    elif (settings.Functional).lower() == 'm062x':
        func1 = '# m062x/'
        func2 = ' int=ultrafine nmr='
        
    else:
        func1 = '# ' + settings.Functional + '/'
        func2 = ' nmr='
            
    if (settings.BasisSet).lower()[:3] == 'pcs':
        basis1 = 'gen'
    else:
        basis1 = settings.BasisSet
            
    CompSettings = func1 + basis1 + func2
    if settings.jJ or settings.jFC:
        CompSettings += '(giao,spinspin,mixed,readatoms)'
    else:
        CompSettings += 'giao'
    
    if settings.Solvent != '':
        CompSettings += ' scrf=(solvent=' + settings.Solvent+')\n'
    else:
        CompSettings += '\n'
        
    f.write(CompSettings)
    f.write('\n'+Gausinp+'\n\n')
    f.write(str(charge) + ' 1\n')

    natom = 0

    for atom in conformer:
        f.write(atoms[natom] + '  ' + atom[1] + '  ' + atom[2] + '  ' +
                atom[3] + '\n')
        natom = natom + 1
    f.write('\n')
    if (settings.BasisSet).lower()[:3] == 'pcs':
        basisfile = file(settings.ScriptDir + '/' + 
                         (settings.BasisSet).lower(), 'r')
        inp = basisfile.readlines()
        basisfile.close()
        for line in inp:
            f.write(line)
        f.write('\n')
    
    if settings.jJ or settings.jFC:
        f.write('atoms=H\n')
    
    f.close()


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


def RunJaguar(folder, JagFiles, settings):
    pass


def IsJagCompleted(f):
    Gfile = open(f, 'r')
    outp = Gfile.readlines()
    Gfile.close()
    if len(outp) < 10:
        return False
    if "Normal termination" in outp[-1]:
        return True
    else:
        return False


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

    print GausNames
    print NTaut
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

            
def main(numDS, MaxEnergy, *allargs):

    args = allargs[:-1]
    ExpNMR = allargs[-1]
    
    RunMacromodel(numDS, *args)    
    
    JagInp = "jaguar run "    
    
    print args    
    
    #Run Jaguar setup script for every diastereomer
    for ds in args:
        SetupInput = 'java Setup Upto ' + ds + ' ' + ds + 'jaginp ' + str(MaxEnergy) + ' 3'
        JagInp = JagInp + ds + 'jaginp*.in '
        print SetupInput
        outp = subprocess.check_output(SetupInput, shell=True)
    
    print JagInp
    
    time.sleep(5)
    
    #Run Jaguar (if there are not output files for all input filed
    #and wait until the last file is completed
    
    inpFiles = glob.glob('*jaginp*.in')
    
    totFiles = len(inpFiles)
    
    compFiles = len(glob.glob('*jaginp*.out'))
    
    if compFiles < totFiles:
        outp = subprocess.check_output(JagInp, shell=True)
    
    print outp    
    
    #Wait while all the .out files have been created and report on progress
    while compFiles < totFiles:
        outpFiles = glob.glob('*jaginp*.out')
        if len(outpFiles) > compFiles:
            compFiles = len(outpFiles)
            if compFiles > 1:
                print "Jaguar calculation for file " + str(compFiles-1) + " out of " + str(totFiles) + " completed."
        time.sleep(30)
        
    lastSeries = glob.glob(args[-1] + '*jaginp*.out')
    numLast = len(lastSeries)
    for f in lastSeries:
        if str(numLast) in f[-7:]:
            lastFile = f
            
    #Wait while the last job has been completed by looking in the file
    lastOutp = open(lastFile, 'r')
    outp = lastOutp.readlines()
    
    while not "completed" in outp[-1]:
        lastOutp.close()
        time.sleep(15)
        lastOutp = open(lastFile, 'r')
        outp = lastOutp.readlines()
    
    lastOutp.close()
        
    print "Jaguar calculation for file " + str(compFiles) + " out of " + str(totFiles) + " completed."
    
    #Run NMRDP4
    import NMRDP4
    
    if (numDS<2):
        print "DP4 requires at least 2 candidate structures!"
    else:
        NMRargs = []
        for ds in args:
            NMRargs.append(ds + 'jaginp')
            NMRargs.append(len(glob.glob(ds + '*jaginp*.out')))
        print numDS
        NMRargs.append(ExpNMR)
        print NMRargs
        NMRDP4.main(numDS, *NMRargs)