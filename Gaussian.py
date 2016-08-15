#!/usr/bin/env python
from __future__ import division
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 15:56:54 2014

@author: ke291

Contains all of the Gaussian specific code for input generation and calculation
execution. Called by PyDP4.py.
"""

import Tinker
import MacroModel
import nmrPredictGaus

import subprocess
import socket
import os
import time
import sys
import glob
import shutil
import math
import pyximport
pyximport.install()
import ConfPrune


def SetupGaussian(MMoutp, Gausinp, numDigits, settings, adjRMSDcutoff):

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


def PM7opt(Gausinp, conformer, atoms, charge, settings):
    
    WriteMopacFile(Gausinp, conformer, atoms, charge)
    #outp = subprocess.check_output(settings.MOPAC + Gausinp + 'm.mop', shell=True)
    outp = subprocess.Popen([settings.MOPAC, Gausinp + 'm.mop'], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    OptConformer = ReadMopacFile(Gausinp+'m.out')
    
    os.rename(Gausinp+'m.out',Gausinp+'.mopout')
    return OptConformer
    

def WriteMopacFile(Gausinp, conformer, atoms, charge):
    
    f = file(Gausinp + 'm.mop', 'w')

    f.write(" AUX LARGE CHARGE=" + str(charge) + " SINGLET\n")
    #f.write(" AUX LARGE CHARGE=" + str(charge) + " SINGLET RECALC=1 DMAX=0.05 RELSCF=0.01\n")
    f.write(Gausinp+'\n\n')
    
    natom = 0

    for atom in conformer:
        f.write(atoms[natom] + '  ' + atom[1] + '  1  ' + atom[2] + '  1  ' +
                atom[3] + '  1 \n')
        natom = natom + 1
    f.write('\n')
    
    f.close()


def ReadMopacFile(filename):

    mofile = open(filename, 'r')
    MopacOutp = mofile.readlines()
    mofile.close()
    
    index = 0
    conformer = []
    
    while not 'COMPUTATION TIME' in MopacOutp[index]:
        index = index + 1
    #Find the NMR shielding calculation section
    while not 'CARTESIAN COORDINATES' in MopacOutp[index]:
        index = index + 1
    
    index = index + 2
    
    #Read shielding constants and labels
    for line in MopacOutp[index:]:
        if len(line) < 3:
            break
        else:
            data = filter(None, line[:-1].split(' '))
            conformer.append([data[0]]+data[2:])

    return conformer


def WriteGausFile(Gausinp, conformer, atoms, charge, settings):

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


def WriteGausFileOpt(Gausinp, conformer, atoms, charge, settings):

    #write the initial DFT geometry optimisation input file first
    f1 = file(Gausinp + 'a.com', 'w')
    if(settings.nProc > 1):
        f1.write('%nprocshared=' + str(settings.nProc) + '\n')
    f1.write('%mem=6000MB\n%chk='+Gausinp + '.chk\n')
    
    if settings.DFTOpt:
        if settings.Solvent != '':
            f1.write('# b3lyp/6-31g(d,p) Opt=(maxcycles=50) scrf=(solvent=' +
                     settings.Solvent+')\n')
        else:
            f1.write('# b3lyp/6-31g(d,p) Opt=(maxcycles=50)\n')
    elif settings.PM6Opt:
        if settings.Solvent != '':
            f1.write('# pm6 Opt=(maxcycles=50) scrf=(solvent=' +
                     settings.Solvent+')\n')
        else:
            f1.write('# pm6 Opt=(maxcycles=50)\n')
    elif settings.HFOpt:
        if settings.Solvent != '':
            f1.write('# rhf/6-31g(d,p) Opt=(maxcycles=50) scrf=(solvent=' +
                     settings.Solvent+')\n')
        else:
            f1.write('# rhf/6-31g(d,p) Opt=(maxcycles=50)\n')
    elif settings.M06Opt:
        if settings.Solvent != '':
            f1.write('# m062x/6-31g(d,p) Opt=(maxcycles=50) scrf=(solvent=' +
                     settings.Solvent+')\n')
        else:
            f1.write('# m062x/6-31g(d,p) Opt=(maxcycles=50)\n')
            

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
    f2.write('%mem=6000MB\n%chk='+Gausinp + '.chk\n')
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


#Still need addition of support for geometry optimisation
def RunOnZiggy(folder, queue, GausFiles, settings):

    print "ziggy GAUSSIAN job submission script\n"

    #Check that folder does not exist, create job folder on ziggy
    outp = subprocess.check_output('ssh ziggy ls', shell=True)
    if folder in outp:
        print "Folder exists on ziggy, choose another folder name."
        return

    outp = subprocess.check_output('ssh ziggy mkdir ' + folder, shell=True)

    #Write the qsub scripts
    for f in GausFiles:
        if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
            and (not settings.M06Opt):
            WriteSubScript(f[:-4], queue, folder, settings)
        else:
            WriteSubScriptOpt(f[:-4], queue, folder, settings)
    print str(len(GausFiles)) + ' .qsub scripts generated'

    #Upload .com files and .qsub files to directory
    print "Uploading files to ziggy..."
    for f in GausFiles:
        if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
            and (not settings.M06Opt):
            outp = subprocess.check_output('scp ' + f +' ziggy:~/' + folder,
                                           shell=True)
        else:
            outp = subprocess.check_output('scp ' + f[:-4] +'a.com ziggy:~/' +
                                           folder, shell=True)
            outp = subprocess.check_output('scp ' + f[:-4] +'b.com ziggy:~/' +
                                           folder, shell=True)
        outp = subprocess.check_output('scp ' + f[:-4] +'.qsub ziggy:~/' +
                                       folder, shell=True)

    print str(len(GausFiles)) + ' .com and .qsub files uploaded to ziggy'

    #Launch the calculations
    for f in GausFiles:
        job = '~/' + folder + '/' + f[:-4]
        outp = subprocess.check_output('ssh ziggy qsub -l mem=6000mb -q ' + queue + ' -o ' +
            job + '.log -e ' + job + '.err ' + job + '.qsub', shell=True)

    print str(len(GausFiles)) + ' jobs submitted to the queue on ziggy'

    outp = subprocess.check_output('ssh ziggy showq', shell=True)
    if settings.user in outp:
        print "Jobs are running on ziggy"

    Jobs2Complete = list(GausFiles)
    n2complete = len(Jobs2Complete)

    #Check and report on the progress of calculations
    while len(Jobs2Complete) > 0:
        JustCompleted = [job for job in Jobs2Complete if
            IsZiggyGComplete(job[:-3] + 'out', folder, settings)]
        Jobs2Complete[:] = [job for job in Jobs2Complete if
             not IsZiggyGComplete(job[:-3] + 'out', folder, settings)]
        if n2complete != len(Jobs2Complete):
            n2complete = len(Jobs2Complete)
            print str(n2complete) + " remaining."

        time.sleep(60)

    #When done, copy the results back
    print "\nCopying the output files back to localhost..."
    print 'ssh ziggy scp /home/' + settings.user + '/' + folder + '/*.out ' +\
        socket.getfqdn() + ':' + os.getcwd()
    outp = subprocess.check_output('ssh ziggy scp /home/' + settings.user +
                                   '/' + folder + '/*.out ' + socket.getfqdn()
                                   + ':' + os.getcwd(), shell=True)

def RunOnDarwin(folder, GausJobs, settings):
    
    print "Darwin GAUSSIAN job submission script\n"
    
    #Check that folder does not exist, create job folder on ziggy
    outp = subprocess.Popen(['ssh', 'darwin', 'ls'], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    print folder
    
    if folder in outp:
        print "Folder exists on Darwin, choose another folder name."
        return

    outp = subprocess.Popen(['ssh', 'darwin', 'mkdir', folder], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    #Write the slurm scripts
    SubFiles = WriteDarwinScripts(GausJobs, settings)
        
    print str(len(SubFiles)) + ' slurm scripts generated'

    #Upload .com files and slurm files to directory
    print "Uploading files to darwin..."
    for f in GausJobs:
        outp = subprocess.Popen(['scp', f,
            'darwin:/home/' + settings.user + '/' + folder],
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    
    for f in SubFiles:
        
        outp = subprocess.Popen(['scp', f,
            'darwin:/home/' + settings.user + '/' + folder], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    print str(len(GausJobs)) + ' .com and ' + str(len(SubFiles)) +\
        ' slurm files uploaded to darwin'
    
    fullfolder = '/home/' + settings.user + '/' + folder
    #Launch the calculations
    for f in SubFiles:
        outp = subprocess.Popen(['ssh', 'darwin', 'cd ' + fullfolder + ';sbatch', f], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        print "\n".split(outp)[-2:]

    print str(len(SubFiles)) + ' jobs submitted to the queue on darwin ' + \
        'containing ' + str(len(GausJobs)) + ' Gaussian jobs'
    
    Jobs2Complete = list(GausJobs)
    n2complete = len(Jobs2Complete)
    
    #Check and report on the progress of calculations
    while len(Jobs2Complete) > 0:
        JobFinished = IsDarwinGComplete(Jobs2Complete, folder, settings)
        
        Jobs2Complete[:] = [job for job in Jobs2Complete if
             not JobFinished[job[:-3] + 'out']]
        if n2complete != len(Jobs2Complete):
            n2complete = len(Jobs2Complete)
            print str(n2complete) + " remaining."

        time.sleep(180)

    #When done, copy the results back
    print "\nCopying the output files back to localhost..."
    print 'scp darwin:' + fullfolder + '/*.out ' + os.getcwd() + '/'
    #outp = subprocess.check_output('ssh darwin scp /home/' + settings.user +
    #                               '/' + folder + '/*.out ' + socket.getfqdn()
    #                               + ':' + os.getcwd(), shell=True)
        
    outp = subprocess.Popen(['scp', 'darwin:' + fullfolder + '/*.out ',
            os.getcwd() + '/'], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    
    print "\nDeleting checkpoint files..."
    print 'ssh darwin rm ' + fullfolder + '/*.chk'
    #outp = subprocess.check_output('ssh darwin rm /home/' + settings.user +
    #                               '/' + folder + '/*.chk', shell=True)
    outp = subprocess.Popen(['ssh', 'darwin', 'rm', fullfolder + '/*.chk'], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    

def WriteDarwinScripts(GausJobs, settings):
    
    SubFiles = []
    NodeSize = settings.DarwinNodeSize
    
    if len(GausJobs) <= NodeSize:

        SubFiles.append(WriteSlurm(GausJobs, settings))
    
    else:
        print "The Gaussian calculations will be submitted as " +\
                    str(math.ceil(len(GausJobs)/NodeSize)) + \
                    " jobs"
        i = 0
        while (i+1)*NodeSize < len(GausJobs):
            PartGausJobs = list(GausJobs[(i*NodeSize):((i+1)*NodeSize)])
            print "Writing script nr " + str(i+1)
            
            SubFiles.append(WriteSlurm(PartGausJobs, settings, str(i+1)))
            
            i += 1
        
        PartGausJobs = list(GausJobs[(i*NodeSize):])
        print "Writing script nr " + str(i+1)
        SubFiles.append(WriteSlurm(PartGausJobs, settings, str(i+1)))
        
    return SubFiles


def WriteSlurm(GausJobs, settings, index=''):
    
    cwd = os.getcwd()
    filename = settings.Title + 'slurm' + index
    
    shutil.copyfile(settings.ScriptDir + '/Defaultslurm',
                    cwd + '/' + filename)
    slurmf = open(filename, 'r+')
    slurm = slurmf.readlines()
    slurm[12] = '#SBATCH -J ' + settings.Title + '\n'
    slurm[18] = '#SBATCH --ntasks=' + str(len(GausJobs)) + '\n'
    slurm[20] = '#SBATCH --time=' + format(settings.TimeLimit,"02") +\
        ':00:00\n'
    
    for f in GausJobs:
        
        slurm.append('srun --exclusive -n 1 $application < ' + f[:-3] + \
            'com > ' + f[:-3] + 'out 2> error &\n')
            
    slurm.append('wait\n')
    slurmf.truncate(0)
    slurmf.seek(0)
    slurmf.writelines(slurm)
    
    return filename


def IsDarwinGComplete(GausJobs, folder, settings):

    path = '/home/' + settings.user + '/' + folder + '/'
    results = {}
    
    for f in GausJobs:
        outp = subprocess.Popen(['ssh', 'darwin', 'cat', path + f[:-3] + 'out'], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        #outp = subprocess.check_output('ssh darwin cat ' + path + f,
        #                                    shell=True)
        if "Normal termination" in outp[-90:]:
            results[f[:-3] + 'out'] = True
        else:
            results[f[:-3] + 'out'] = False
    
    return results


def WriteSubScript(GausJob, queue, ZiggyJobFolder, settings):

    if not (os.path.exists(GausJob+'.com')):
        print "The input file " + GausJob + ".com does not exist. Exiting..."
        return

    #Create the submission script
    QSub = open(GausJob + ".qsub", 'w')

    #Choose the queue
    QSub.write('#PBS -q ' + queue + '\n')
    if settings.nProc >1:
        QSub.write('#PBS -l nodes=1:ppn=' + str(settings.nProc) + '\n#\n')
    else:
        QSub.write('#PBS -l nodes=1:ppn=1\n#\n')

    #define input files and output files
    QSub.write('file=' + GausJob + '\n\n')
    QSub.write('inpfile=${file}.com\noutfile=${file}.out\n')

    #define cwd and scratch folder and ask the machine
    #to make it before running the job
    QSub.write('HERE=/home/' + settings.user +'/' + ZiggyJobFolder + '\n')
    QSub.write('LSCRATCH=/scratch/' + settings.user + '/' +
               GausJob + '\n')
    QSub.write('mkdir ${LSCRATCH}\n')

    #Setup GAUSSIAN environment variables
    if settings.nProc >1:
        QSub.write('set OMP_NUM_THREADS=12\n')
    else:
        QSub.write('set OMP_NUM_THREADS=1\n')
    QSub.write('export GAUSS_EXEDIR=/usr/local/shared/gaussian/em64t/09-D01/g09\n')
    QSub.write('export g09root=/usr/local/shared/gaussian/em64t/09-D01\n')
    QSub.write('export PATH=/usr/local/shared/gaussian/em64t/09-D01/g09:$PATH\n')
    QSub.write('export GAUSS_SCRDIR=$LSCRATCH\n')

    #copy the input file to scratch
    QSub.write('cp ${HERE}/${inpfile}  $LSCRATCH\ncd $LSCRATCH\n')

    #write useful info to the job output file (not the gaussian)
    QSub.write('echo "Starting job $PBS_JOBID"\necho\n')
    QSub.write('echo "PBS assigned me this node:"\ncat $PBS_NODEFILE\necho\n')

    QSub.write('ln -s $HERE/$outfile $LSCRATCH/$outfile\n')
    QSub.write('$GAUSS_EXEDIR/g09 < $inpfile > $outfile\n')

    #Cleanup
    QSub.write('rm -rf ${LSCRATCH}/\n')
    QSub.write('qstat -f $PBS_JOBID\n')

    QSub.close()


#Function to write ziggy script when dft optimisation is used
def WriteSubScriptOpt(GausJob, queue, ZiggyJobFolder, settings):

    if not (os.path.exists(GausJob+'a.com')):
        print "The input file " + GausJob + "a.com does not exist. Exiting..."
        return
    if not (os.path.exists(GausJob+'b.com')):
        print "The input file " + GausJob + "b.com does not exist. Exiting..."
        return

    #Create the submission script
    QSub = open(GausJob + ".qsub", 'w')

    #Choose the queue
    QSub.write('#PBS -q ' + queue + '\n')
    if settings.nProc >1:
        QSub.write('#PBS -l nodes=1:ppn=' + str(settings.nProc) + '\n#\n')
    else:
        QSub.write('#PBS -l nodes=1:ppn=1\n#\n')

    #define input files and output files
    QSub.write('file=' + GausJob + '\n\n')
    QSub.write('inpfile1=${file}a.com\ninpfile2=${file}b.com\n')
    QSub.write('outfile1=${file}temp.out\noutfile2=${file}.out\n')

    #define cwd and scratch folder and ask the machine
    #to make it before running the job
    QSub.write('HERE=/home/' + settings.user +'/' + ZiggyJobFolder + '\n')
    QSub.write('LSCRATCH=/scratch/' + settings.user + '/' + GausJob + ZiggyJobFolder +'\n')
    QSub.write('mkdir ${LSCRATCH}\n')

    #Setup GAUSSIAN environment variables
    if settings.nProc >1:
        QSub.write('set OMP_NUM_THREADS=12\n')
    else:
        QSub.write('set OMP_NUM_THREADS=1\n')
    QSub.write('export GAUSS_EXEDIR=/usr/local/shared/gaussian/em64t/09-D01/g09\n')
    QSub.write('export g09root=/usr/local/shared/gaussian/em64t/09-D01\n')
    QSub.write('export PATH=/usr/local/shared/gaussian/em64t/09-D01/g09:$PATH\n')
    QSub.write('export GAUSS_SCRDIR=$LSCRATCH\n')

    #copy the input files to scratch
    QSub.write('cp ${HERE}/${inpfile1}  $LSCRATCH\n')
    QSub.write('cp ${HERE}/${inpfile2}  $LSCRATCH\ncd $LSCRATCH\n')

    #write useful info to the job output file (not the gaussian)
    QSub.write('echo "Starting job $PBS_JOBID"\necho\n')
    QSub.write('echo "PBS assigned me this node:"\ncat $PBS_NODEFILE\necho\n')
    
    QSub.write('ln -s $HERE/$outfile1 $LSCRATCH/$outfile1\n')
    QSub.write('ln -s $HERE/$outfile2 $LSCRATCH/$outfile2\n')
    QSub.write('$GAUSS_EXEDIR/g09 < $inpfile1 > $outfile1\n')
    QSub.write('$GAUSS_EXEDIR/g09 < $inpfile2 > $outfile2\n')

    #Cleanup
    QSub.write('rm -rf ${LSCRATCH}/\n')
    QSub.write('qstat -f $PBS_JOBID\n')

    QSub.close()


def IsZiggyGComplete(f, folder, settings):

    path = '/home/' + settings.user + '/' + folder + '/'
    try:
        outp1 = subprocess.check_output('ssh ziggy ls ' + path, shell=True)
    except subprocess.CalledProcessError, e:
        print "ssh ziggy ls failed: " + str(e.output)
        return False
    if f in outp1:
        try:
            outp2 = subprocess.check_output('ssh ziggy cat ' + path + f,
                                            shell=True)
        except subprocess.CalledProcessError, e:
            print "ssh ziggy cat failed: " + str(e.output)
            return False
        if "Normal termination" in outp2[-90:]:
            return True
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
        WriteGausFileOpt(f[:-8], coords,atoms,charge,settings)
        print f[:-8] + '* deleted and new .com files written'
    f = file('Reoptimized.log', 'w')
    f.write('\n'.join([x[:-8] for x in GoutpFiles]))
    f.close()


def ReadGeometry(GOutpFile):

    gausfile = open(GOutpFile, 'r')
    GOutp = gausfile.readlines()

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
            data = filter(None, line[:-1].split(','))
            atoms.append(data[0][-1])
            coords.append(data[1:])
            
    line = GOutp[chindex].split('Charge =  ')
    line = line[1].split(' Multiplicity = ')
    charge = int(line[0])
    print charge
    gausfile.close()

    return atoms, coords, charge


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


def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
