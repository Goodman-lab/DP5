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
import shutil
import math

import Gaussian

SetupNMRCalcs = Gaussian.SetupNMRCalcs

SetupECalcs = Gaussian.SetupECalcs

SetupOptCalcs = Gaussian.SetupOptCalcs

ReadShieldings = Gaussian.ReadShieldings

ReadDFTEnergies = Gaussian.ReadDFTEnergies

IsGausCompleted = Gaussian.IsGausCompleted

def RunNMRCalcs(Isomers, settings):
    print('\nRunning Gaussian NMR calculations on Ziggy...')

    # Run Gaussian jobs on Ziggy cluster in folder named after date and time
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
    print('\nRunning Gaussian DFT energy calculations on Ziggy...')

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


def RunCalcs(GausJobs, settings):

    MaxCon = settings.MaxConcurrentJobs

    if len(GausJobs) < MaxCon:
        if len(GausJobs) > 0:
            RunBatchOnZiggy(0, settings.queue, GausJobs, settings)
    else:
        if len(GausJobs) > 0:
            print("The DFT calculations will be done in " + \
                  str(math.ceil(len(GausJobs) / MaxCon)) + " batches")
            i = 0
            while (i + 1) * MaxCon < len(GausJobs):
                print("Starting batch nr " + str(i + 1))
                RunBatchOnZiggy(str(i + 1), settings.queue, GausJobs[(i * MaxCon):((i + 1) * MaxCon)], settings)
                i += 1
            print("Starting batch nr " + str(i + 1))
            RunBatchOnZiggy(str(i + 1), settings.queue, GausJobs[(i * MaxCon):], settings)

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
