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
import socket
import os
import time
import glob
import shutil
import math


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


def SetupGaussian(Gausinp, settings):

    if settings.charge is not None:
        charge = settings.charge

    for num in range(0, len(pruned)):
        filename = Gausinp+str(num+1).zfill(3)
        if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
            and (not settings.M06Opt):
            WriteGausFile(filename, pruned[num], atoms, charge, settings)
        else:
            if not(os.path.exists(filename + 'temp.out')):
                WriteGausFileOpt(filename, pruned[num], atoms, charge, settings)
            else:
                print('Reusing geometry from ' + filename + 'temp.out')
                tempatoms, tempcoords, tempcharge = ReadTempGeometry(filename + 'temp.out')
                if (len(tempatoms) > 0) and (len(tempcoords) > 0):
                    WriteGausFileOpt(filename, tempcoords, tempatoms, tempcharge, settings)
                else:
                    print(filename + 'temp.out is an invalid Gaussian output file. Falling back to MMFF geometry.')
                    WriteGausFileOpt(filename, pruned[num], atoms, charge, settings)

    print(str(len(pruned)) + " .com files written")


def WriteGausFile(Gausinp, conformer, atoms, charge, settings):

    f = open(Gausinp + '.com', 'w')
    if(settings.nProc > 1):
        f.write('%nprocshared=' + str(settings.nProc) + '\n')
    if settings.DFT == 'g':
        f.write('%mem=2000MB\n%chk='+Gausinp + '.chk\n')
    else:
        f.write('%mem=6000MB\n%chk='+Gausinp + '.chk\n')
    
    if (settings.nFunctional).lower() == 'wp04':
        func1 = '# blyp/'
        func2 = ' iop(3/76=1000001189,3/77=0961409999,3/78=0000109999)' + \
            ' nmr='
        
    elif (settings.nFunctional).lower() == 'm062x':
        func1 = '# m062x/'
        func2 = ' int=ultrafine nmr='
        
    else:
        func1 = '# ' + settings.nFunctional + '/'
        func2 = ' nmr='
            
    if (settings.nBasisSet).lower()[:3] == 'pcs':
        basis1 = 'gen'
    else:
        basis1 = settings.nBasisSet
            
    CompSettings = func1 + basis1 + func2
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
        f.write(atoms[natom] + '  ' + atom[0] + '  ' + atom[1] + '  ' +
                atom[2] + '\n')
        natom = natom + 1
    f.write('\n')
    if (settings.nBasisSet).lower()[:3] == 'pcs':
        basisfile = file(settings.ScriptDir + '/' + 
                         (settings.nBasisSet).lower(), 'r')
        inp = basisfile.readlines()
        basisfile.close()
        for line in inp:
            f.write(line)
        f.write('\n')
    
    f.close()


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


#Still need addition of support for geometry optimisation
def RunOnZiggy(findex, queue, GausFiles, settings):

    if findex == 0:
        folder = settings.StartTime + settings.Title
    else:
        folder = settings.StartTime + findex + settings.Title

    scrfolder = settings.StartTime + settings.Title

    print("ziggy GAUSSIAN job submission script\n")

    #Check that folder does not exist, create job folder on ziggy
    outp = subprocess.check_output('ssh ziggy ls', shell=True)
    if folder in outp:
        print("Folder exists on ziggy, choose another folder name.")
        return

    outp = subprocess.check_output('ssh ziggy mkdir ' + folder, shell=True)

    print("Results folder: " + folder)

    #Write the qsub scripts
    for f in GausFiles:
        if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
            and (not settings.M06Opt):
            WriteSubScript(f[:-4], queue, folder, settings)
        else:
            WriteSubScriptOpt(f[:-4], queue, folder, settings)
    print(str(len(GausFiles)) + ' slurm scripts generated')

    #Upload .com files and .qsub files to directory
    print("Uploading files to ziggy...")
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
        outp = subprocess.check_output('scp ' + f[:-4] +'slurm ziggy:~/' +
                                       folder, shell=True)

    print(str(len(GausFiles)) + ' .com and slurm files uploaded to ziggy')

    #Launch the calculations
    for f in GausFiles:
        job = '~/' + folder + '/' + f[:-4]
        outp = subprocess.check_output('ssh ziggy sbatch ' + job + 'slurm', shell=True)

    print(str(len(GausFiles)) + ' jobs submitted to the queue on ziggy')

    outp = subprocess.check_output('ssh ziggy qstat', shell=True)
    if settings.user in outp:
        print("Jobs are running on ziggy")

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
            print(str(n2complete) + " remaining.")

        time.sleep(60)

    #When done, copy the results back
    print("\nCopying the output files back to localhost...")
    print('ssh ziggy scp /home/' + settings.user + '/' + folder + '/*.out ' +\
        socket.getfqdn() + ':' + os.getcwd())
    outp = subprocess.check_output('ssh ziggy scp /home/' + settings.user +
                                   '/' + folder + '/*.out ' + socket.getfqdn()
                                   + ':' + os.getcwd(), shell=True)

def RunOnDarwin(findex, GausJobs, settings):

    if findex == 0:
        folder = settings.StartTime + settings.Title
    else:
        folder = settings.StartTime + findex + settings.Title

    scrfolder = settings.StartTime + settings.Title

    print("Darwin GAUSSIAN job submission script\n")
    
    #Check that results folder does not exist, create job folder on darwin
    outp = subprocess.Popen(['ssh', 'darwin', 'ls'], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    print("Results folder: " + folder)
    
    if folder in outp:
        print("Results folder exists on Darwin, choose another folder name.")
        quit()

    outp = subprocess.Popen(['ssh', 'darwin', 'mkdir', folder], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    # Check that scratch directory does not exist, create job folder on darwin
    outp = subprocess.Popen(['ssh', 'darwin', 'ls ' + settings.DarwinScrDir], \
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    print("Scratch directory: " + settings.DarwinScrDir + scrfolder)

    if folder in outp:
        print("Scratch folder exists on Darwin, choose another folder name.")
        quit()

    outp = subprocess.Popen(['ssh', 'darwin', 'mkdir', settings.DarwinScrDir + scrfolder], \
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    #Write the slurm scripts
    SubFiles = WriteDarwinScripts(GausJobs, settings, scrfolder)
        
    print(str(len(SubFiles)) + ' slurm scripts generated')

    #Upload .com files and slurm files to directory
    print("Uploading files to darwin...")
    for f in GausJobs:
        if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
            and (not settings.M06Opt):
            outp = subprocess.Popen(['scp', f,
            'darwin:/home/' + settings.user + '/' + folder],
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        
        else:
            outp = subprocess.Popen(['scp', f[:-4] + 'a.com',
                'darwin:/home/' + settings.user + '/' + folder],
                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            outp = subprocess.Popen(['scp', f[:-4] + 'b.com',
                'darwin:/home/' + settings.user + '/' + folder],
                stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    
    for f in SubFiles:
        
        outp = subprocess.Popen(['scp', f,
            'darwin:/home/' + settings.user + '/' + folder], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    print(str(len(GausJobs)) + ' .com and ' + str(len(SubFiles)) +\
        ' slurm files uploaded to darwin')
    
    fullfolder = '/home/' + settings.user + '/' + folder
    #Launch the calculations
    for f in SubFiles:
        outp = subprocess.Popen(['ssh', 'darwin', 'cd ' + fullfolder + ';sbatch', f], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        print(outp.split('\n')[-2])

    print(str(len(SubFiles)) + ' jobs submitted to the queue on darwin ' + \
        'containing ' + str(len(GausJobs)) + ' Gaussian jobs')
    
    Jobs2Complete = list(GausJobs)
    n2complete = len(Jobs2Complete)
    
    #Check and report on the progress of calculations
    while len(Jobs2Complete) > 0:
        JobFinished = IsDarwinGComplete(Jobs2Complete, folder, settings)
        
        Jobs2Complete[:] = [job for job in Jobs2Complete if
             not JobFinished[job[:-3] + 'out']]
        if n2complete != len(Jobs2Complete):
            n2complete = len(Jobs2Complete)
            print(str(n2complete) + " remaining.")

        time.sleep(180)

    #When done, copy the results back
    print("\nCopying the output files back to localhost...")
    print('scp darwin:' + fullfolder + '/*.out ' + os.getcwd() + '/')
    #outp = subprocess.check_output('ssh darwin scp /home/' + settings.user +
    #                               '/' + folder + '/*.out ' + socket.getfqdn()
    #                               + ':' + os.getcwd(), shell=True)
        
    outp = subprocess.Popen(['scp', 'darwin:' + fullfolder + '/*.out ',
            os.getcwd() + '/'], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    fullscrfolder = settings.DarwinScrDir + scrfolder
    print("\nDeleting scratch folder...")
    print('ssh darwin rm -r ' + fullscrfolder)
    #outp = subprocess.check_output('ssh darwin rm /home/' + settings.user +
    #                               '/' + folder + '/*.chk', shell=True)
    outp = subprocess.Popen(['ssh', 'darwin', 'rm -r', fullscrfolder], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    

def WriteDarwinScripts(GausJobs, settings, scrfolder):

    SubFiles = []
    NodeSize = settings.DarwinNodeSize
    AdjNodeSize = int(math.floor(settings.DarwinNodeSize/settings.nProc))

    if len(GausJobs) <= AdjNodeSize:

        SubFiles.append(WriteSlurm(GausJobs, settings, scrfolder))
    
    else:
        print("The Gaussian calculations will be submitted as " +\
                    str(math.ceil(len(GausJobs)/AdjNodeSize)) + \
                    " jobs")
        i = 0
        while (i+1)*AdjNodeSize < len(GausJobs):
            PartGausJobs = list(GausJobs[(i*AdjNodeSize):((i+1)*AdjNodeSize)])
            print("Writing script nr " + str(i+1))
            
            SubFiles.append(WriteSlurm(PartGausJobs, settings, scrfolder, str(i+1)))
            
            i += 1
        
        PartGausJobs = list(GausJobs[(i*AdjNodeSize):])
        print("Writing script nr " + str(i+1))
        SubFiles.append(WriteSlurm(PartGausJobs, settings, scrfolder, str(i+1)))
        
    return SubFiles


def WriteSlurm(GausJobs, settings, scrfolder, index=''):
    
    cwd = os.getcwd()
    filename = settings.Title + 'slurm' + index
    
    shutil.copyfile(settings.ScriptDir + '/Defaultslurm',
                    cwd + '/' + filename)
    slurmf = open(filename, 'r+')
    slurm = slurmf.readlines()
    slurm[12] = '#SBATCH -J ' + settings.Title + '\n'
    slurm[19] = '#SBATCH --ntasks=' + str(len(GausJobs)*settings.nProc) + '\n'
    slurm[21] = '#SBATCH --time=' + format(settings.TimeLimit,"02") +\
        ':00:00\n'

    slurm[61] = 'export GAUSS_SCRDIR=' + settings.DarwinScrDir + scrfolder + '\n'
    
    if (not settings.DFTOpt) and (not settings.PM6Opt) and (not settings.HFOpt)\
        and (not settings.M06Opt):
            
        for f in GausJobs:
            slurm.append('srun --exclusive -n1 $application < ' + f[:-3] + \
                'com > ' + f[:-3] + 'out 2> error &\n')
            #slurm.append('$application < ' + f[:-3] + \
            #             'com > ' + f[:-3] + 'out 2> error &\n')
        slurm.append('wait\n')
    else:
        for f in GausJobs:
            slurm.append('(srun --exclusive -n1 -c' + str(settings.nProc) + ' $application < ' + f[:-4] + \
                'a.com > ' + f[:-4] + 'temp.out 2> error;')
            slurm.append('srun --exclusive -n1 -c' + str(settings.nProc) + ' $application < ' + f[:-4] + \
                         'b.com > ' + f[:-4] + '.out 2> error) &\n')
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
        print("The input file " + GausJob + ".com does not exist. Exiting...")
        return

    #Create the submission script
    QSub = open(GausJob + "slurm", 'w')

    #Choose the queue
    QSub.write('#!/bin/bash\n\n')
    QSub.write('#SBATCH -p SWAN\n')
    if settings.nProc >1:
        QSub.write('#SBATCH --nodes=1\n#SBATCH --cpus-per-task=' + str(settings.nProc) + '\n')
    else:
        QSub.write('#SBATCH --nodes=1\n#SBATCH --cpus-per-task=1\n')
    QSub.write('#SBATCH --time=' + format(settings.TimeLimit,"02") +':00:00\n\n')

    #define input files and output files
    QSub.write('file=' + GausJob + '\n\n')
    QSub.write('inpfile=${file}.com\noutfile=${file}.out\n')

    #define cwd and scratch folder and ask the machine
    #to make it before running the job
    QSub.write('HERE=/home/' + settings.user +'/' + ZiggyJobFolder + '\n')
    QSub.write('SCRATCH=/scratch/' + settings.user + '/' +
               GausJob + '\n')
    QSub.write('mkdir ${SCRATCH}\n')

    #Setup GAUSSIAN environment variables
    QSub.write('set OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n')
    
    QSub.write('export GAUSS_EXEDIR=/usr/local/shared/gaussian/em64t/09-D01/g09\n')
    QSub.write('export g09root=/usr/local/shared/gaussian/em64t/09-D01\n')
    QSub.write('export PATH=/usr/local/shared/gaussian/em64t/09-D01/g09:$PATH\n')
    QSub.write('export GAUSS_SCRDIR=$SCRATCH\n')
    QSub.write('exe=$GAUSS_EXEDIR/g09\n')
    #copy the input file to scratch
    QSub.write('cp ${HERE}/${inpfile}  $SCRATCH\ncd $SCRATCH\n')

    #write useful info to the job output file (not the gaussian)
    QSub.write('echo "Starting job $SLURM_JOBID"\necho\n')
    QSub.write('echo "SLURM assigned me this node:"\nsrun hostname\necho\n')

    QSub.write('ln -s $HERE/$outfile $SCRATCH/$outfile\n')
    QSub.write('srun $exe > $outfile < $inpfile\n')

    #Cleanup
    QSub.write('rm -rf ${SCRATCH}/\n')
    
    QSub.close()


#Function to write ziggy script when dft optimisation is used
def WriteSubScriptOpt(GausJob, queue, ZiggyJobFolder, settings):

    if not (os.path.exists(GausJob+'a.com')):
        print("The input file " + GausJob + "a.com does not exist. Exiting...")
        return
    if not (os.path.exists(GausJob+'b.com')):
        print("The input file " + GausJob + "b.com does not exist. Exiting...")
        return

    #Create the submission script
    QSub = open(GausJob + "slurm", 'w')

    #Choose the queue
    QSub.write('#!/bin/bash\n\n')
    QSub.write('#SBATCH -p SWAN\n')
    if settings.nProc >1:
        QSub.write('#SBATCH --nodes=1\n#SBATCH --cpus-per-task=' + str(settings.nProc) + '\n')
    else:
        QSub.write('#SBATCH --nodes=1\n#SBATCH --cpus-per-task=1\n')
    QSub.write('#SBATCH --time=' + format(settings.TimeLimit,"02") +':00:00\n\n')

    #define input files and output files
    QSub.write('file=' + GausJob + '\n\n')
    QSub.write('inpfile1=${file}a.com\ninpfile2=${file}b.com\n')
    QSub.write('outfile1=${file}temp.out\noutfile2=${file}.out\n')

    #define cwd and scratch folder and ask the machine
    #to make it before running the job
    QSub.write('HERE=/home/' + settings.user +'/' + ZiggyJobFolder + '\n')
    QSub.write('SCRATCH=/scratch/' + settings.user + '/' +
               GausJob + '\n')
    QSub.write('mkdir ${SCRATCH}\n')

    #Setup GAUSSIAN environment variables
    QSub.write('set OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n')
    
    QSub.write('export GAUSS_EXEDIR=/usr/local/shared/gaussian/em64t/09-D01/g09\n')
    QSub.write('export g09root=/usr/local/shared/gaussian/em64t/09-D01\n')
    QSub.write('export PATH=/usr/local/shared/gaussian/em64t/09-D01/g09:$PATH\n')
    QSub.write('export GAUSS_SCRDIR=$SCRATCH\n')
    QSub.write('exe=$GAUSS_EXEDIR/g09\n')

    #copy the input files to scratch
    QSub.write('cp ${HERE}/${inpfile1}  $SCRATCH\n')
    QSub.write('cp ${HERE}/${inpfile2}  $SCRATCH\ncd $SCRATCH\n')

    #write useful info to the job output file (not the gaussian)
    QSub.write('echo "Starting job $SLURM_JOBID"\necho\n')
    QSub.write('echo "SLURM assigned me this node:"\nsrun hostname\necho\n')
    
    QSub.write('ln -s $HERE/$outfile1 $SCRATCH/$outfile1\n')
    QSub.write('ln -s $HERE/$outfile2 $SCRATCH/$outfile2\n')
    QSub.write('srun $exe < $inpfile1 > $outfile1\n')
    QSub.write('wait\n')
    QSub.write('srun $exe < $inpfile2 > $outfile2\n')

    #Cleanup
    QSub.write('rm -rf ${SCRATCH}/\n')
    
    QSub.close()


def IsZiggyGComplete(f, folder, settings):

    path = '/home/' + settings.user + '/' + folder + '/'
    try:
        outp1 = subprocess.check_output('ssh ziggy ls ' + path, shell=True)
    except subprocess.CalledProcessError as e:
        print("ssh ziggy ls failed: " + str(e.output))
        return False
    if f in outp1:
        try:
            outp2 = subprocess.check_output('ssh ziggy cat ' + path + f,
                                            shell=True)
        except subprocess.CalledProcessError as e:
            print("ssh ziggy cat failed: " + str(e.output))
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