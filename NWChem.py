# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 15:56:54 2014

@author: ke291

Contains all of the NWChem specific code for input generation and calculation
execution. Called by PyDP4.py.
"""

import Tinker
import MacroModel
import nmrPredictNWChem

import glob
import os
import subprocess
import time
import socket
"""
if os.name == 'nt':
    import pyximport
    pyximport.install()
    import ConfPrune
else:
    import pyximport
    pyximport.install()
    import ConfPrune
"""

def SetupNWChem(MMoutp, NWCheminp, numDigits, settings, adjRMSDcutoff):

    #Reads conformer geometry, energies and atom labels from Tinker output
    #(atoms, conformers) = ReadConformers(MMoutp, MaxEnergy)

    if settings.MMTinker:
        #Reads conformer geometry, energies and atom labels from Tinker output
        (atoms, conformers, charge) = Tinker.ReadTinker(MMoutp, settings)
    else:
        (atoms, conformers, charge) = MacroModel.ReadMacromodel(MMoutp,
                                                                settings)
    if settings.charge is not None:
        charge = settings.charge
    #Prune similar conformations, if the number exceeds the limit
    if len(conformers) > settings.PerStructConfLimit:
        pruned = ConfPrune.RMSDPrune(conformers, atoms, adjRMSDcutoff)
    else:
        pruned = conformers

    print str(len(conformers) - len(pruned)) +\
        " or " + "{:.1f}".format(100 * (len(conformers) - len(pruned)) /
        len(conformers)) + "% of conformations have been pruned based on " + \
        str(adjRMSDcutoff) + " angstrom cutoff"

    for num in range(0, len(pruned)):
        filename = NWCheminp+str(num+1).zfill(3)
        WriteNWChemFile(filename, pruned[num], atoms, charge, settings)

    print str(len(pruned)) + " .nw files written"


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


def WriteNWChemFile(NWCheminp, conformer, atoms, charge, settings):

    f = file(NWCheminp + '.nw', 'w')
    f.write('memory stack 1500 mb heap 1500 mb global 3000 mb\n')
    if settings.DFT == 'w':
        f.write('scratch_dir /scratch/' + settings.user + '/' + NWCheminp + '\n')
    f.write('echo\n\nstart molecule\n\ntitle "'+NWCheminp+'"\n')
    f.write('echo\n\nstart\n\n')
    
    if settings.charge is not None:
        f.write('charge ' + str(settings.charge) + '\n\n')
    else:
        f.write('charge ' + str(charge) + '\n\n')
    
    f.write('geometry units angstroms print xyz autosym\n')

    natom = 0
    for atom in conformer:
        f.write('  ' + atoms[natom] + ' ' + atom[1] + ' ' + atom[2] + ' ' +
                atom[3] + '\n')
        natom = natom + 1
    
    basis = settings.BasisSet
    if basis.lower() == '6-31g(d,p)':
        basis = '6-31g**'
    elif basis.lower() == '6-311g(d)':
        basis = '6-311g*'
    
    f.write('end\n\nbasis\n  * library ' + basis + '\nend\n\n')
    if settings.Solvent != "":
        f.write('cosmo\n  do_cosmo_smd true\n  solvent ' + settings.Solvent + '\n')
        f.write('end\n\n')
    if settings.DFTOpt or settings.HFOpt:
        f.write('driver\n  maxiter ' + str(settings.MaxDFTOptCycles)+ '\nend\n\n')
        if settings.DFTOpt:
            f.write('dft\n  xc b3lyp\n  mult 1\nend\n\n')
            f.write('task dft optimize\n\n')
        if settings.M06Opt:
            f.write('dft\n  xc m06-2x\n  mult 1\nend\n\n')
            f.write('task dft optimize\n\n')
        if settings.HFOpt:
            f.write('task scf optimize\n\n')
    if (settings.Functional).lower() == 'b3lyp':
        f.write('dft\n  xc b3lyp\n  mult 1\nend\n\n')
    elif (settings.Functional).lower() == 'm062x' or\
        (settings.Functional).lower() == 'm06-2x':
        f.write('dft\n  xc m06-2x\n  mult 1\nend\n\n')
    elif (settings.Functional).lower() == 'mpw1pw91':
        f.write('dft\n  xc mpw91 0.75 HFexch 0.25 perdew91\n  mult 1\nend\n\n')
    else:
        f.write('dft\n  xc ' + settings.Functional + '\n  mult 1\nend\n\n')
    f.write('task dft energy\n\n')
    f.write('property\n  shielding\nend\n')
    f.write('task dft property\n')
    f.close()


def GetFiles2Run(inpfiles, settings):
    #Get the names of all relevant input files
    NinpFiles = []
    for filename in inpfiles:
        NinpFiles = NinpFiles + glob.glob(filename + 'nwinp???.nw')

    Files2Run = []

    #for every input file check that there is a completed output file,
    #delete the incomplete outputs and add the inputs to be done to Files2Run
    for filename in NinpFiles:
        if not os.path.exists(filename[:-3]+'.nwo'):
            Files2Run.append(filename)
        else:
            if IsNWChemCompleted(filename[:-3] + '.nwo'):
                continue
            else:
                os.remove(filename[:-3] + '.nwo')
                Files2Run.append(filename)

    return Files2Run


def RunNWChem(inpfiles, settings):

    NCompleted = 0
    NWChemPrefix = "nwchem "

    for f in inpfiles:
        print NWChemPrefix + f + ' > ' + f[:-2] + 'nwo'
        outp = subprocess.check_output(NWChemPrefix + f + ' > ' + f[:-2] +
                                       'nwo', shell=True)
        NCompleted += 1
        print "NWChem job " + str(NCompleted) + " of " + str(len(inpfiles)) + \
            " completed."


def RunNMRPredict(numDS, *args):

    NWNames = []
    NTaut = []

    for val in range(0, numDS):
        NTaut.append(args[val*2])
        NWNames.append(args[val*2+1])

    RelEs = []
    populations = []
    BoltzmannShieldings = []
    SigConfs = []

    print NWNames
    print NTaut
    #This loop runs nmrPredict for each diastereomer and collects
    #the outputs    
    for isomer in NWNames:

        NWFiles = glob.glob(isomer + 'nwinp*.nwo')
        for f in range(0, len(NWFiles)):
            NWFiles[f] = NWFiles[f][:-4]

        #Runs nmrPredictNWChem Name001, ... and collects output
        (x, y, labels, z, SCs) = nmrPredictNWChem.main(*NWFiles)
        RelEs.append(x)
        populations.append(y)
        BoltzmannShieldings.append(z)
        SigConfs.append(SCs)

    return (RelEs, populations, labels, BoltzmannShieldings, SigConfs, NTaut)


def IsNWChemCompleted(f):
    Nfile = open(f, 'r')
    outp = Nfile.readlines()
    Nfile.close()
    outp = "".join(outp)
    if "AUTHORS" in outp:
        return True
    else:
        return False


def RunOnZiggy(folder, queue, NWFiles, settings):

    print "ziggy NWChem job submission script\n"

    #Check that folder does not exist, create job folder on ziggy
    outp = subprocess.check_output('ssh ziggy ls', shell=True)
    if folder in outp:
        print "Folder exists on ziggy, choose another folder name."
        return

    outp = subprocess.check_output('ssh ziggy mkdir ' + folder, shell=True)
    #Write the qsub scripts
    for f in NWFiles:
        WriteSubScript(f[:-3], queue, folder, settings)
    print str(len(NWFiles)) + ' .qsub scripts generated'

    #Upload .com files and .qsub files to directory
    print "Uploading files to ziggy..."
    for f in NWFiles:
        outp = subprocess.check_output('scp ' + f +' ziggy:~/' + folder,
                                       shell=True)
        outp = subprocess.check_output('scp ' + f[:-3] +'.qsub ziggy:~/' +
                                       folder, shell=True)

    print str(len(NWFiles)) + ' .nw and .qsub files uploaded to ziggy'

    #Launch the calculations
    for f in NWFiles:
        job = '~/' + folder + '/' + f[:-3]
        outp = subprocess.check_output('ssh ziggy qsub -q ' + queue + ' -o '
            + job + '.log -e ' + job + '.err -l nodes=1:ppn=1:ivybridge ' +
            job + '.qsub', shell=True)
        time.sleep(3)

    print str(len(NWFiles)) + ' jobs submitted to the queue on ziggy'

    outp = subprocess.check_output('ssh ziggy showq', shell=True)
    if settings.user in outp:
        print "Jobs are running on ziggy"

    Jobs2Complete = list(NWFiles)
    n2complete = len(Jobs2Complete)

    #Check and report on the progress of calculations
    while len(Jobs2Complete) > 0:
        JustCompleted = [job for job in Jobs2Complete if
            IsZiggyGComplete(job[:-2] + 'nwo', folder, settings)]
        Jobs2Complete[:] = [job for job in Jobs2Complete if
             not IsZiggyGComplete(job[:-2] + 'nwo', folder, settings)]
        if n2complete != len(Jobs2Complete):
            n2complete = len(Jobs2Complete)
            print str(n2complete) + " remaining."

        time.sleep(60)

    #When done, copy the results back
    print "\nCopying the output files back to localhost..."
    print 'ssh ziggy scp /home/' + settings.user + '/' + folder + '/*.nwo ' +\
        socket.getfqdn() + ':' + os.getcwd()
    outp = subprocess.check_output('ssh ziggy scp /home/' + settings.user + '/'
                                   + folder + '/*.nwo ' + socket.getfqdn() + ':'
                                   + os.getcwd(), shell=True)


def RunOnMedivir(NWFiles, settings):

    print "Medivir NWChem job submission script\n"

    #Write the qsub scripts
    for f in NWFiles:
        WriteMedivirSubScript(f[:-3], settings)
    print str(len(NWFiles)) + ' .qsub scripts generated'

    #Launch the calculations
    for f in NWFiles:
        job = f[:-3]
        outp = subprocess.check_output('qsub ' + job + '.qsub', shell=True)
        time.sleep(3)

    print str(len(NWFiles)) + ' jobs submitted to the queue on ziggy'

    outp = subprocess.check_output('qstat', shell=True)
    if 'nwinp' in outp:
        print "Jobs are running on the cluster"

    Jobs2Complete = list(NWFiles)
    n2complete = len(Jobs2Complete)

    #Check and report on the progress of calculations
    while len(Jobs2Complete) > 0:
        JustCompleted = [job for job in Jobs2Complete if
            IsMedivirComplete(job[:-2] + 'nwo', settings)]
        Jobs2Complete[:] = [job for job in Jobs2Complete if
             not IsMedivirComplete(job[:-2] + 'nwo', settings)]
        if n2complete != len(Jobs2Complete):
            n2complete = len(Jobs2Complete)
            print str(n2complete) + " remaining."

        time.sleep(60)

    print "Calculation on the cluster done.\n"
                                   
                                   
def WriteSubScript(NWJob, queue, ZiggyJobFolder, settings):

    if not (os.path.exists(NWJob+'.nw')):
        print "The input file " + NWJob + ".nw does not exist. Exiting..."
        return

    #Create the submission script
    QSub = open(NWJob + ".qsub", 'w')

    #Choose the queue
    QSub.write('#PBS -q ' + queue + '\n#PBS -l nodes=1:ppn=1\n#\n')

    #define input files and output files
    QSub.write('file=' + NWJob + '\n\n')
    QSub.write('inpfile=${file}.nw\noutfile=${file}.nwo\n')

    #define cwd and scratch folder and ask the machine
    #to make it before running the job
    QSub.write('HERE=/home/' + settings.user + '/' + ZiggyJobFolder + '\n')
    QSub.write('SCRATCH=/sharedscratch/' + settings.user + '/' + NWJob + '\n')
    QSub.write('LSCRATCH=/scratch/' + settings.user + '/' + NWJob + '\n')
    QSub.write('mkdir ${SCRATCH}\n')
    QSub.write('mkdir ${LSCRATCH}\n')

    #load relevant modules
    QSub.write('set OMP_NUM_THREADS=1\n')
    QSub.write('module load anaconda\n')
    QSub.write('module load gcc/4.8.3\n')
    QSub.write('module load mpi/openmpi/gnu/1.8.1\n')
    QSub.write('module load nwchem\n')

    #copy the input file to scratch
    QSub.write('cp ${HERE}/${inpfile}  $SCRATCH\ncd $SCRATCH\n')

    #write useful info to the job output file (not the gaussian)
    QSub.write('echo "Starting job $PBS_JOBID"\necho\n')
    QSub.write('echo "PBS assigned me this node:"\ncat $PBS_NODEFILE\necho\n')

    QSub.write('ln -s $HERE/$outfile $SCRATCH/$outfile\n')
    QSub.write('nwchem $inpfile > $outfile\n')

    #Cleanup
    QSub.write('rm -rf ${SCRATCH}/\n')
    QSub.write('rm -rf ${LSCRATCH}/\n')
    QSub.write('qstat -f $PBS_JOBID\n')

    QSub.close()


def WriteMedivirSubScript(NWJob, settings):

    if not (os.path.exists(NWJob+'.nw')):
        print "The input file " + NWJob + ".nw does not exist. Exiting..."
        return

    #Create the submission script
    QSub = open(NWJob + ".qsub", 'w')

    QSub.write('#PBS -S /bin/tcsh\n')
    QSub.write('#PBS -l nodes=1:ppn=24\n#\n')
    QSub.write('#PBS -N ' + NWJob + '\n')
    QSub.write('#PBS -j oe\n')

    QSub.write('mkdir /scr/' + NWJob + '\n')
    QSub.write('cd /scr/' + NWJob + '\n')

    #load relevant modules
    QSub.write('module load nwchem-6.5\n')

    QSub.write('mpiexec nwchem ' + os.getcwd() + '/' + NWJob + '.nw > '
	 + os.getcwd() + '/' + NWJob + '.nwo\n')

    QSub.write('cd ../\n')
    QSub.write('rm -r /scr/' + NWJob + '\n')

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
        if "AUTHORS" in outp2:
            return True
    return False


def IsMedivirComplete(f, settings):

    try:
        outp1 = subprocess.check_output('ls ', shell=True)
    except subprocess.CalledProcessError, e:
        print "ls failed: " + str(e.output)
        return False
    if f in outp1:
        try:
            outp2 = subprocess.check_output('cat ' + f, shell=True)
        except subprocess.CalledProcessError, e:
            print "cat failed: " + str(e.output)
            return False
        if "AUTHORS" in outp2:
            return True
    return False
