#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 15:56:54 2019

@author: ke291

Contains all of the specific code for running Gaussian jobs on Darwin.
A lot of code is reused from Gaussian.py. Called by PyDP4.py.
"""

import subprocess
import os
import time
import shutil
import math

import Gaussian

MaxConcurrentJobs = 160

SetupNMRCalcs = Gaussian.SetupNMRCalcs

SetupECalcs = Gaussian.SetupECalcs

SetupOptCalcs = Gaussian.SetupOptCalcs

ReadShieldings = Gaussian.ReadShieldings

ReadEnergies = Gaussian.ReadEnergies

ReadGeometries = Gaussian.ReadGeometries

IsGausCompleted = Gaussian.IsGausCompleted

Converged = Gaussian.Converged

def RunNMRCalcs(Isomers, settings):
    print('\nRunning Gaussian NMR calculations on Darwin...')

    # Run Gaussian jobs on CSD3 cluster in folder named after date and time
    # Split in batches, if needed

    jobdir = os.getcwd()
    os.chdir('nmr')

    GausJobs = []

    for iso in Isomers:
        GausJobs.extend([x for x in iso.NMRInputFiles if (x[:-4] + '.out') not in iso.NMROutputFiles])

    Completed = RunCalcs(GausJobs, settings)

    for iso in Isomers:
        iso.NMROutputFiles.extend([x[:-4] + '.out' for x in iso.NMRInputFiles if (x[:-4] + '.out') in Completed])

    os.chdir(jobdir)

    return Isomers


def RunECalcs(Isomers, settings):
    print('\nRunning Gaussian DFT energy calculations on Darwin...')

    jobdir = os.getcwd()
    os.chdir('e')

    GausJobs = []

    for iso in Isomers:
        GausJobs.extend([x for x in iso.EInputFiles if (x[:-4] + '.out') not in iso.EOutputFiles])

    Completed = RunCalcs(GausJobs, settings)

    for iso in Isomers:
        iso.EOutputFiles.extend([x[:-4] + '.out' for x in iso.EInputFiles if (x[:-4] + '.out') in Completed])

    os.chdir(jobdir)

    return Isomers


def RunOptCalcs(Isomers, settings):
    print('\nRunning Gaussian DFT geometry optimizations on Darwin...')

    jobdir = os.getcwd()
    os.chdir('opt')

    GausJobs = []

    for iso in Isomers:
        GausJobs.extend([x for x in iso.OptInputFiles if (x[:-4] + '.out') not in iso.OptOutputFiles])

    Completed = RunCalcs(GausJobs, settings)

    for iso in Isomers:
        iso.OptOutputFiles.extend([x[:-4] + '.out' for x in iso.OptInputFiles if (x[:-4] + '.out') in Completed])

    os.chdir(jobdir)

    return Isomers


def RunCalcs(GausJobs, settings):

    MaxCon = MaxConcurrentJobs

    if len(GausJobs) < MaxCon:
        if len(GausJobs) > 0:
            RunBatchOnDarwin(0, GausJobs, settings)
    else:
        if len(GausJobs) > 0:
            print("The DFT calculations will be done in " + \
                  str(math.ceil(len(GausJobs) / MaxCon)) + " batches")
            i = 0
            while (i + 1) * MaxCon < len(GausJobs):
                print("Starting batch nr " + str(i + 1))
                RunBatchOnDarwin(str(i + 1), GausJobs[(i * MaxCon):((i + 1) * MaxCon)], settings)
                i += 1
            print("Starting batch nr " + str(i + 1))
            RunBatchOnDarwin(str(i + 1), GausJobs[(i * MaxCon):], settings)

    NCompleted = 0
    Completed = []

    for f in GausJobs:
        if IsGausCompleted(f[:-4] + '.out'):
            Completed.append(f[:-4] + '.out')
            NCompleted += 1

    print(str(NCompleted) + "Gaussian jobs of " + str(len(GausJobs)) + \
        " completed successfully.")

    return Completed


def RunBatchOnDarwin(findex, GausJobs, settings):

    if findex == 0:
        folder = settings.StartTime + settings.Title
    else:
        folder = settings.StartTime + findex + settings.Title

    scrfolder = settings.StartTime + settings.Title

    #Check that results folder does not exist, create job folder on darwin
    outp = subprocess.Popen(['ssh', 'darwin', 'ls'], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    print("Results folder: " + folder)
    
    if folder in outp.decode():
        print("Results folder exists on Darwin, choose another folder name.")
        quit()

    outp = subprocess.Popen(['ssh', 'darwin', 'mkdir', folder], \
      stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    # Check that scratch directory does not exist, create job folder on darwin
    outp = subprocess.Popen(['ssh', 'darwin', 'ls ' + settings.DarwinScrDir], \
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    print("Scratch directory: " + settings.DarwinScrDir + scrfolder)

    if folder in outp.decode():
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
        outp = subprocess.Popen(['scp', f,
        'darwin:/home/' + settings.user + '/' + folder],
        stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        
    for f in SubFiles:
        
        outp = subprocess.Popen(['scp', f,
            'darwin:/home/' + settings.user + '/' + folder], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    print(str(len(GausJobs)) + ' .com and ' + str(len(SubFiles)) +\
        ' slurm files uploaded to darwin')
    
    fullfolder = '/home/' + settings.user + '/' + folder
    JobIDs = []

    #Launch the calculations
    for f in SubFiles:
        outp = subprocess.Popen(['ssh', 'darwin', 'cd ' + fullfolder + ';sbatch', f], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        status = outp.decode().split('\n')[-2]
        print(status)
        JobIDs.append(status.split('job ')[1])

    print(str(len(SubFiles)) + ' jobs submitted to the queue on darwin ' + \
        'containing ' + str(len(GausJobs)) + ' Gaussian jobs')

    time.sleep(60)

    OldQRes = CheckDarwinQueue(JobIDs, settings)

    while OldQRes[0] < 0:
        OldQRes = CheckDarwinQueue(JobIDs, settings)
        time.sleep(60)

    print('Pending: ' + str(OldQRes[0]) + ', Running: ' + str(OldQRes[1]) + ', Not in queue: ' + str(OldQRes[2]))

    Jobs2Complete = list(GausJobs)
    n2complete = len(Jobs2Complete)
    
    #Check and report on the progress of calculations
    while len(Jobs2Complete) > 0:
        JobFinished = IsDarwinGComplete(Jobs2Complete, folder, settings)
        
        Jobs2Complete[:] = [job for job in Jobs2Complete if
             not JobFinished[job[:-3] + 'out']]
        if n2complete != len(Jobs2Complete):
            n2complete = len(Jobs2Complete)
            print(str(n2complete) + " Gaussian jobs remaining.")

        QRes = CheckDarwinQueue(JobIDs, settings)
        if QRes != OldQRes:
            if QRes[0] < 0:
                QRes = OldQRes
            else:
                OldQRes = QRes
                print('Darwin queue:')
                print('Pending: ' + str(OldQRes[0]) + ', Running: ' + str(OldQRes[1]) + ', Not in queue: ' + str(OldQRes[2]))


        if (QRes[2] == len(JobIDs)) and (QRes[0] >= 0):
            #check each gaussian file to ascertain the status of individual gaus jobs
            print('No jobs left in Darwin queue')
            break

        time.sleep(180)

    #When done, copy the results back
    print("\nCopying the output files back to localhost...")
    print('scp darwin:' + fullfolder + '/*.out ' + os.getcwd() + '/')

    outp = subprocess.Popen(['scp', 'darwin:' + fullfolder + '/*.out',
            os.getcwd() + '/'], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    print("\nDeleting checkpoint files...")
    print('ssh darwin rm ' + fullfolder + '/*.chk')
    outp = subprocess.Popen(['ssh', 'darwin', 'rm', fullfolder + '/*.chk'], \
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]

    fullscrfolder = settings.DarwinScrDir + scrfolder
    print("\nDeleting scratch folder...")
    print('ssh darwin rm -r ' + fullscrfolder)

    outp = subprocess.Popen(['ssh', 'darwin', 'rm -r', fullscrfolder], \
            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    

def WriteDarwinScripts(GausJobs, settings, scrfolder):

    SubFiles = []
    AdjNodeSize = int(math.floor(settings.DarwinNodeSize/settings.nProc))

    # If the jobs exactly fill the node, just write the submission script
    if len(GausJobs) == AdjNodeSize:
        SubFiles.append(WriteSlurm(GausJobs, settings, scrfolder))

    # If the jobs do not fill the node, increase the processor count to fill it
    elif len(GausJobs) < AdjNodeSize:
        NewNProc = int(math.floor(settings.DarwinNodeSize / len(GausJobs)))
        SubFiles.append(WriteSlurm(GausJobs, settings, scrfolder, nProc=NewNProc))
        print("Jobs don't fill the Darwin node, nproc increased to " + str(NewNProc))
        for j, GausJob in enumerate(GausJobs):
            line = '%nprocshared=' + str(NewNProc) + '\n'
            ReplaceLine(GausJob, 0, line)

    # If the jobs more than fill the node, submit as several jobs
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

        # if the last few jobs do not fill the node, increase the processor count to fill it
        if len(PartGausJobs) < AdjNodeSize:
            NewNProc = int(math.floor(settings.DarwinNodeSize / len(PartGausJobs)))
            print("Jobs don't fill the last Darwin node, nproc increased to " + str(NewNProc))
            print("Writing script nr " + str(i + 1))
            SubFiles.append(WriteSlurm(PartGausJobs, settings, scrfolder, str(i+1), nProc=NewNProc))
            for j, GausJob in enumerate(PartGausJobs):
                line = '%nprocshared=' + str(NewNProc) + '\n'
                ReplaceLine(GausJob, 0, line)
        else:
            print("Writing script nr " + str(i + 1))
            SubFiles.append(WriteSlurm(PartGausJobs, settings, scrfolder, str(i+1)))

    return SubFiles


def WriteSlurm(GausJobs, settings, scrfolder, index='', nProc = -1):

    if nProc == -1:
        nProc = settings.nProc

    cwd = os.getcwd()
    filename = settings.Title + 'slurm' + index
    
    shutil.copyfile(settings.ScriptDir + '/Defaultslurm',
                    cwd + '/' + filename)
    slurmf = open(filename, 'r+')
    slurm = slurmf.readlines()
    slurm[12] = '#SBATCH -J ' + settings.Title + '\n'
    slurm[14] = '#SBATCH -A ' + settings.project + '\n'
    slurm[19] = '#SBATCH --ntasks=' + str(len(GausJobs)*nProc) + '\n'
    slurm[21] = '#SBATCH --time=' + format(settings.TimeLimit,"02") +\
        ':00:00\n'

    slurm[61] = 'export GAUSS_SCRDIR=' + settings.DarwinScrDir + scrfolder + '\n'
    
    for f in GausJobs:
        slurm.append('srun --exclusive -n1 -c' + str(nProc) + ' $application < ' + f[:-3] + \
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

        if "Normal termination" in outp[-90:].decode():
            results[f[:-3] + 'out'] = True
        else:
            results[f[:-3] + 'out'] = False
    
    return results


def CheckDarwinQueue(JobIDs, settings):

    outp = subprocess.Popen(['ssh', 'darwin', 'squeue', '-u ' + settings.user], \
                            stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    outp = outp.decode().split('\n')
    QStart = -1000
    for i, line in enumerate(outp):
        if 'JOBID' in line:
            QStart = i+1
            break

    if QStart < 0:
        return -100, -100, -100

    QueueReport = outp[QStart:-1]
    JobStats = []

    for job in JobIDs:
        status = ''
        for i, line in enumerate(QueueReport):
            if job in line:
                status = list(filter(None, line.split(' ')))[4]
        JobStats.append(status)

    Pending = JobStats.count('PD')
    Running = JobStats.count('R')
    NotInQueue = JobStats.count('')

    return Pending, Running, NotInQueue


def ReplaceLine(File, LineN, Line):
    gausf = open(File, 'r+')
    gauslines = gausf.readlines()
    gauslines[LineN] = Line
    gausf.truncate(0)
    gausf.seek(0)
    gausf.writelines(gauslines)
    gausf.close()
