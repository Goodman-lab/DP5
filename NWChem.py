# -*- coding: utf-8 -*-
"""
Rewritten on Wed Jan 15 2020

@author: ke291

Contains all of the NWChem specific code for input generation and calculation
execution. Called by PyDP4.py.
"""

import glob
import os
import subprocess
import shutil

def SetupNMRCalcs(Isomers, settings):
    jobdir = os.getcwd()

    if not os.path.exists('nmr'):
        os.mkdir('nmr')
    os.chdir('nmr')

    for iso in Isomers:
        if iso.ExtCharge > -10:
            charge = iso.ExtCharge
        else:
            charge = iso.MMCharge

        if iso.DFTConformers == []:
            conformers = iso.Conformers
        else:
            conformers = iso.DFTConformers

        for num in range(0, len(conformers)):
            filename = iso.BaseName + 'nwinp' + str(num + 1).zfill(3)

            if os.path.exists(filename + '.nwo'):
                if IsNWChemCompleted(filename + '.nwo'):
                    iso.NMROutputFiles.append(filename + '.nwo')
                    continue
                else:
                    os.remove(filename + '.nwo')

            WriteNWChemFile(filename, conformers[num], iso.Atoms, charge, settings, 'nmr')
            iso.NMRInputFiles.append(filename + '.nw')

    os.chdir(jobdir)

    return Isomers


def SetupECalcs(Isomers, settings):

    jobdir = os.getcwd()

    if not os.path.exists('e'):
        os.mkdir('e')
    os.chdir('e')

    for iso in Isomers:
        if iso.ExtCharge > -10:
            charge = iso.ExtCharge
        else:
            charge = iso.MMCharge

        if iso.DFTConformers == []:
            conformers = iso.Conformers
        else:
            conformers = iso.DFTConformers

        for num in range(0, len(conformers)):
            filename = iso.BaseName + 'nwinp' + str(num + 1).zfill(3)

            if os.path.exists(filename + '.nwo'):
                if IsNWChemCompleted(filename + '.nwo'):
                    iso.EOutputFiles.append(filename + '.nwo')
                    continue
                else:
                    os.remove(filename + '.nwo')

            WriteNWChemFile(filename, conformers[num], iso.Atoms, charge, settings, 'e')
            iso.EInputFiles.append(filename + '.nw')

    os.chdir(jobdir)

    return Isomers


def SetupOptCalcs(Isomers, settings):
    jobdir = os.getcwd()

    if not os.path.exists('opt'):
        os.mkdir('opt')
    os.chdir('opt')

    for iso in Isomers:
        if iso.ExtCharge > -10:
            charge = iso.ExtCharge
        else:
            charge = iso.MMCharge

        if iso.DFTConformers == []:
            conformers = iso.Conformers
        else:
            conformers = iso.DFTConformers

        for num in range(0, len(conformers)):
            filename = iso.BaseName + 'nwinp' + str(num + 1).zfill(3)

            if os.path.exists(filename + '.nwo'):
                if IsNWChemCompleted(filename + '.nwo'):
                    iso.OptOutputFiles.append(filename + '.nwo')
                    continue
                else:
                    os.remove(filename + '.nwo')

            WriteNWChemFile(filename, conformers[num], iso.Atoms, charge, settings, 'opt')
            iso.OptInputFiles.append(filename + '.nw')

    os.chdir(jobdir)

    return Isomers


def Converged(Isomers):

    jobdir = os.getcwd()

    if not os.path.exists('opt'):
        os.chdir(jobdir)
        return False

    # insert code for convergence testing here

    os.chdir('opt')

    os.chdir(jobdir)
    return True


def RunNMRCalcs(Isomers, settings):

    print('\nRunning NWChem DFT NMR calculations locally...')

    jobdir = os.getcwd()
    os.chdir('nmr')

    NWJobs = []

    for iso in Isomers:
        print(iso.NMRInputFiles)
        NWJobs.extend([x for x in iso.NMRInputFiles if (x[:-3] + '.nwo') not in iso.NMROutputFiles])

    Completed = RunCalcs(NWJobs, settings)

    for iso in Isomers:
        iso.NMROutputFiles.extend([x[:-3] + '.nwo' for x in iso.NMRInputFiles if (x[:-3] + '.nwo') in Completed])

    os.chdir(jobdir)

    return Isomers


def GetPrerunNMRCalcs(Isomers):

    print('\nLooking for prerun NWChem DFT NMR files...')

    jobdir = os.getcwd()
    os.chdir('nmr')
    """
    for iso in Isomers:
        iso.NMRInputFiles = glob.glob(iso.BaseName + 'ginp*com')
        iso.NMROutputFiles.extend([x[:-4] + '.out' for x in iso.NMRInputFiles if IsGausCompleted(x[:-4] + '.out')])
    """
    print('NMR calc files:')
    print(', '.join([', '.join(x.NMROutputFiles) for x in Isomers]))

    os.chdir(jobdir)

    return Isomers


def RunECalcs(Isomers, settings):

    print('\nRunning NWChem DFT energy calculations locally...')

    jobdir = os.getcwd()
    os.chdir('e')

    NWJobs = []
    for iso in Isomers:
        print(iso.EInputFiles)
        NWJobs.extend([x for x in iso.EInputFiles if (x[:-3] + '.nwo') not in iso.EOutputFiles])

    print('To run: ' + str(NWJobs))
    Completed = RunCalcs(NWJobs, settings)

    for iso in Isomers:
        iso.EOutputFiles.extend([x[:-3] + '.nwo' for x in iso.EInputFiles if (x[:-3] + '.nwo') in Completed])

    os.chdir(jobdir)

    return Isomers


def GetPrerunECalcs(Isomers):

    print('\nLooking for prerun NWChem DFT energy calculation files...')

    jobdir = os.getcwd()
    os.chdir('e')

    """
    for iso in Isomers:
        iso.EInputFiles = glob.glob(iso.BaseName + 'ginp*com')
        iso.EOutputFiles.extend([x[:-4] + '.out' for x in iso.EInputFiles if IsGausCompleted(x[:-4] + '.out')])
    """
    print('Energy files:')
    print(', '.join([', '.join(x.EOutputFiles) for x in Isomers]))

    os.chdir(jobdir)

    return Isomers


def RunOptCalcs(Isomers, settings):

    print('\nRunning NWChem DFT geometry optimizations locally...')

    jobdir = os.getcwd()
    os.chdir('opt')

    NWJobs = []

    for iso in Isomers:
        print(iso.OptInputFiles)
        NWJobs.extend([x for x in iso.OptInputFiles if (x[:-3] + '.nwo') not in iso.OptOutputFiles])

    Completed = RunCalcs(NWJobs, settings)

    for iso in Isomers:
        iso.OptOutputFiles.extend([x[:-3] + '.nwo' for x in iso.OptInputFiles if (x[:-3] + '.nwo') in Completed])

    os.chdir(jobdir)

    return Isomers


def GetPrerunOptCalcs(Isomers):

    print('\nLooking for prerun NWChem DFT optimization files...')

    jobdir = os.getcwd()
    os.chdir('opt')
    """
    for iso in Isomers:
        iso.OptInputFiles = glob.glob(iso.BaseName + 'ginp*com')
        iso.OptOutputFiles.extend([x[:-4] + '.out' for x in iso.OptInputFiles if IsGausCompleted(x[:-4] + '.out')])
    """
    print('Opt files:')
    print(', '.join([', '.join(x.OptOutputFiles) for x in Isomers]))

    os.chdir(jobdir)

    return Isomers


def RunCalcs(NWJobs, settings):

    NCompleted = 0
    Completed = []
    NWChemPrefix = settings.NWChemPath

    if shutil.which(NWChemPrefix) is None:
        print('NWChem.py, RunCalcs:\n  Could not find NWChem executable at ' + NWChemPrefix)
        quit()

    for f in NWJobs:
        print(NWChemPrefix + ' ' + f + ' > ' + f[:-2] + 'nwo')
        outp = subprocess.check_output(NWChemPrefix + ' ' + f + ' > ' + f[:-2] +
                                       'nwo', shell=True)
        NCompleted += 1
        print("NWChem job " + str(NCompleted) + " of " + str(len(NWJobs)) + \
            " completed.")
        if IsNWChemCompleted(f[:-3] + '.nwo'):
            Completed.append(f[:-3] + '.nwo')
            print("NWChem job " + str(NCompleted) + " of " + str(len(NWJobs)) + \
                  " completed.")
        else:
            print("NWChem job terminated with an error. Continuing.")

    if NCompleted > 0:
        print(str(NCompleted) + " NWChem jobs completed successfully.")
    elif len(NWJobs) == 0:
        print("There were no jobs to run.")

    return Completed


def WriteNWChemFile(NWinp, conformer, atoms, charge, settings, type):

    f = open(NWinp + '.nw', 'w')
    f.write('memory stack 1500 mb heap 1500 mb global 3000 mb\n')
    if settings.DFT == 'w':
        f.write('scratch_dir /scratch/' + settings.user + '/' + NWinp + '\n')
    f.write('echo\n\nstart molecule\n\ntitle "' + NWinp + '"\n')
    f.write('echo\n\nstart\n\n')

    if settings.charge is not None:
        f.write('charge ' + str(settings.charge) + '\n\n')
    else:
        f.write('charge ' + str(charge) + '\n\n')

    f.write('geometry units angstroms print xyz autosym\n')

    natom = 0
    for atom in conformer:
        f.write('  ' + atoms[natom] + ' ' + atom[0] + ' ' + atom[1] + ' ' +
                atom[2] + '\n')
        natom = natom + 1

    basis = settings.nBasisSet
    if basis.lower() == '6-31g(d,p)':
        basis = '6-31g**'
    elif basis.lower() == '6-311g(d)':
        basis = '6-311g*'

    f.write('end\n\nbasis\n  * library ' + basis + '\nend\n\n')
    if settings.Solvent != "":
        GausSolvents = ['chloroform', 'dimethylsulfoxide', 'benzene', 'methanol', 'pyridine', 'acetone']
        NWSolvents = ['chcl3', 'dmso', 'benzene', 'methanol', 'pyridine', 'acetone']

        if settings.Solvent in GausSolvents:
            solvent = NWSolvents[GausSolvents.index(settings.Solvent)]
        else:
            solvent = settings.Solvent

        f.write('cosmo\n  do_cosmo_smd true\n  solvent ' + solvent + '\n')
        f.write('end\n\n')

    if type == 'nmr':
        f.write(NMRSuffix(settings))
    elif type == 'e':
        f.write(ESuffix(settings))
    elif type == 'opt':
        f.write(OptSuffix(settings))

    f.close()


def NMRSuffix(settings):
    suffix = ''
    if (settings.nFunctional).lower() == 'b3lyp':
        suffix += 'dft\n  xc b3lyp\n  mult 1\nend\n\n'
    elif (settings.nFunctional).lower() == 'm062x' or \
            (settings.nFunctional).lower() == 'm06-2x':
        suffix +=  'dft\n  xc m06-2x\n  mult 1\nend\n\n'
    elif (settings.nFunctional).lower() == 'mpw1pw91':
        suffix += 'dft\n  xc mpw91 0.75 HFexch 0.25 perdew91\n  mult 1\nend\n\n'
    else:
        suffix += 'dft\n  xc ' + settings.nFunctional + '\n  mult 1\nend\n\n'
    suffix += 'task dft energy\n\nproperty\n  shielding\nend\ntask dft property\n'

    return suffix


def ESuffix(settings):
    suffix = ''
    if (settings.nFunctional).lower() == 'b3lyp':
        suffix += 'dft\n  xc b3lyp\n  mult 1\nend\n\n'
    elif (settings.nFunctional).lower() == 'm062x' or \
            (settings.nFunctional).lower() == 'm06-2x':
        suffix +=  'dft\n  xc m06-2x\n  mult 1\nend\n\n'
    elif (settings.nFunctional).lower() == 'mpw1pw91':
        suffix += 'dft\n  xc mpw91 0.75 HFexch 0.25 perdew91\n  mult 1\nend\n\n'
    else:
        suffix += 'dft\n  xc ' + settings.nFunctional + '\n  mult 1\nend\n\n'
    suffix += 'task dft energy\n'

    return suffix


def OptSuffix(settings):
    suffix = 'driver\n  maxiter ' + str(settings.MaxDFTOptCycles) + '\nend\n\n'
    if (settings.oFunctional).lower() == 'b3lyp':
        suffix += 'dft\n  xc b3lyp\n  mult 1\nend\n\n'
        suffix += 'task dft optimize\n\n'
    elif (settings.oFunctional).lower() == 'm062x' or (settings.oFunctional).lower() == 'm06-2x':
        suffix += 'dft\n  xc m06-2x\n  mult 1\nend\n\n'
        suffix += 'task dft optimize\n\n'

    return suffix


def IsNWChemCompleted(f):
    Nfile = open(f, 'r')
    outp = Nfile.readlines()
    Nfile.close()
    outp = "".join(outp)
    if "AUTHORS" in outp:
        return True
    else:
        return False


def IsNWChemConverged(f):
    Nfile = open(f, 'r')
    outp = Nfile.readlines()
    Nfile.close()
    outp = "".join(outp)
    if 'Optimization converged' in outp:
        return True
    else:
        return False


#Read energy from e, if not present, then o, if not present, then nmr
def ReadEnergies(Isomers, settings):
    jobdir = os.getcwd()

    if 'e' in settings.Workflow:
        os.chdir('e')
    elif 'o' in settings.Workflow:
        os.chdir('opt')
    else:
        os.chdir('nmr')

    for i, iso in enumerate(Isomers):

        if 'e' in settings.Workflow:
            NWOutpFiles = iso.EOutputFiles
        elif 'o' in settings.Workflow:
            NWOutpFiles = iso.OptOutputFiles
        else:
            NWOutpFiles = iso.NMROutputFiles

        DFTEnergies = []
        for NWOutpFile in NWOutpFiles:
            nwfile = open(NWOutpFile, 'r')
            NWOutp = nwfile.readlines()
            nwfile.close()

            for line in NWOutp:
                if 'Total DFT energy' in line:
                    start = line.index('Total')
                    energy = float(line[start + 19:])

            DFTEnergies.append(energy)

        Isomers[i].DFTEnergies = DFTEnergies

    os.chdir(jobdir)
    return Isomers


def ReadShieldings(Isomers):
    jobdir = os.getcwd()
    os.chdir('nmr')

    for iso in Isomers:

        if len(iso.NMROutputFiles) < 1:
            print("NWChem.py, ReadShieldings: No NMR DFT output" +
                  " files found, NMR data could not be read. Quitting.")
            quit()

        for NWOutpFile in iso.NMROutputFiles:
            nwfile = open(NWOutpFile, 'r')
            NWOutp = nwfile.readlines()
            nwfile.close()

            index = 0
            shieldings = []
            labels = []

            # Find the NMR shielding calculation section
            while not 'Chemical Shielding' in NWOutp[index]:
                index = index + 1

            # Read shielding constants and labels
            for line in NWOutp[index:]:
                if 'isotropic' in line:
                    start = line.index('isotropic')
                    shieldings.append(float(line[start + 13:]))
                if 'Atom' in line:
                    start = line.index('Atom')
                    labels.append(line[start + 12] + line[start + 7:start + 10].strip())

            print(NWOutpFile,len(shieldings))

            iso.ConformerShieldings.append(shieldings)

        iso.ShieldingLabels = labels

    os.chdir(jobdir)

    return Isomers


def ReadGeometry(NWOutpFile):

    nwfile = open(NWOutpFile, 'r')
    NWOutp = nwfile.readlines()
    nwfile.close()

    atoms = []
    coords = []
    gindex = -1

    # Find the last geometry section
    for index in range(len(NWOutp)):
        if ('Geometry "geometry" -> "geometry"' in NWOutp[index]):
            gindex = index + 7

    if gindex < 0:
        print('Error: No geometry found in file ' + NWOutpFile)
        quit()

    # Read geometry
    for line in NWOutp[gindex:]:
        if len(line) < 2:
            break
        else:
            data = [_f for _f in line[:-1].split(' ') if _f]
            atoms.append(data[1])
            coords.append(data[3:])

    return atoms, coords


def ReadGeometries(Isomers, settings):

    jobdir = os.getcwd()
    if ('o' in settings.Workflow):
        os.chdir('opt')

        for iso in Isomers:

            iso.DFTConformers = [[] for x in iso.OptOutputFiles]

            if len(iso.OptOutputFiles) < 1:
                print("NWChem.py, ReadGeometries: No geometry optimisation output" +
                      " files found, geometries could not be read. Quitting.")
                quit()

            for num, NWOutpFile in enumerate(iso.OptOutputFiles):

                atoms, coords = ReadGeometry(NWOutpFile)

                iso.DFTConformers[num] = coords

            iso.Atoms = atoms
    else:
        os.chdir('nmr')

        for iso in Isomers:

            iso.DFTConformers = [[] for x in iso.NMROutputFiles]

            if len(iso.NMROutputFiles) < 1:
                print("NWChem.py, ReadGeometries: No geometry optimisation output" +
                      " files found, geometries could not be read. Quitting.")
                quit()

            for num, NWOutpFile in enumerate(iso.NMROutputFiles):
                atoms, coords = ReadGeometry(NWOutpFile)

                iso.DFTConformers[num] = coords

            iso.Atoms = atoms
    #return atoms, coords, charge
    os.chdir(jobdir)
    return Isomers

