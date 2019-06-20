# -*- coding: utf-8 -*-
"""
Created on Mon Jan  12 14:42:47 2015

@author: ke291

Takes care of all the NMR description interpretation, equivalent atom
averaging, Boltzmann averaging and DP4 input preparation and running DP4.py. Called by PyDP4.py

FUNCTIONS AFTER REWRITE:
Calculation of NMR shifts based on TMS reference
Equivalent atom averaging
NMR description parsing
NMR raw data interpretation top level organization
"""

import re
import os
import math

gasConstant = 8.3145
temperature = 298.15
hartreeEnergy = 2625.499629554010

# Data struture for loading and keeping all of experimental NMR data in one place.
class NMRData:
    def __init__(self, InputPath):
        self.InputPath = InputPath  # Initial structure input file
        self.Type = 'desc'          # desc or fid, depending on whether the description or raw data used
        self.Atoms = []             # Element labels
        self.Cshifts = []           # Experimental C NMR shifts
        self.Clabels = []           # Experimental C NMR labels, if any
        self.Hshifts = []           # Experimental H NMR shifts
        self.Hlabels = []           # Experimental H NMR labels, if any
        self.Equivalents = []       # Atoms assumed to be NMR equivalent in computational data
        self.Omits = []

        if not os.path.exists(self.InputPath):
            print('NMR data path does not exist, quitting...')
            quit()
        if os.path.isfile(self.InputPath):
            self.Type = 'desc'
            self.ExpNMRFromDesc()
        else:
            self.Type = 'fid'
            print('FID data reading and analysis not implemented yet, quitting...')
            quit()


    def ExpNMRFromDesc(self):

        print('Loading NMR data from ' + self.InputPath)

        # Reads the experimental NMR data from the file
        ExpNMR_file = open(self.InputPath, 'r')
        Cexp = ExpNMR_file.readline()
        ExpNMR_file.readline()
        Hexp = ExpNMR_file.readline()

        # Check if exp NMR file contains info about equivalent atoms and read it
        # into an array
        # Also reads a list of atoms to omit from analysis

        equivalents = []
        omits = []

        ExpNMR_file.readline()
        for line in ExpNMR_file:
            if not 'OMIT' in line and len(line) > 1:
                equivalents.append(line[:-1].split(','))
            elif 'OMIT' in line:
                omits.extend(line[5:-1].split(','))

        ExpNMR_file.close()

        self.Clabels, self.Cshifts = self.ParseExp(Cexp)
        self.Hlabels, self.Hshifts = self.ParseExp(Hexp)
        self.Equivalents = equivalents
        self.Omits = omits


    def ParseExp(self, exp):
        # Replace all 'or' and 'OR' with ',', remove all spaces and 'any'
        texp = re.sub(r"or|OR", ',', exp, flags=re.DOTALL)
        texp = re.sub(r" |any", '', texp, flags=re.DOTALL)

        # Get all assignments, split mulitassignments
        expLabels = re.findall(r"(?<=\().*?(?=\)|;)", texp, flags=re.DOTALL)
        expLabels = [x.split(',') for x in expLabels]

        # Remove assignments and get shifts
        ShiftData = (re.sub(r"\(.*?\)", "", exp, flags=re.DOTALL)).split(',')
        expShifts = [float(x) for x in ShiftData]

        return expLabels, expShifts


def CalcBoltzmannWeightedShieldings(Isomers):

    energies = []

    for i, iso in enumerate(Isomers):

        # Calculate rel. energies in kJ/mol
        minE = min(iso.DFTEnergies)

        relEs = []

        for e in iso.DFTEnergies:
            relEs.append((e - minE) * hartreeEnergy)

        Isomers[i].Energies = relEs

        populations = []

        # Calculate Boltzmann populations
        for e in relEs:
            populations.append(math.exp(-e * 1000 / (gasConstant * temperature)))

        q = sum(populations)

        for p in range(0, len(populations)):
            populations[p] = populations[p] / q

        Isomers[i].Populations = populations

        print(Isomers[i].Populations)

        # Calculate Boltzmann weighed shielding constants
        # by summing the shifts multiplied by the isomers population
        BoltzmannShieldings = []

        for atom in range(len(iso.Atoms)):
            shielding = 0
            for population, shieldings in zip(iso.Populations, iso.ConformerShieldings):
                shielding = shielding + shieldings[atom] * population
            BoltzmannShieldings.append(shielding)

        Isomers[i].BoltzmannShieldings = BoltzmannShieldings

    return Isomers


def GetTMSConstants(settings):
    TMSfile = open(settings.ScriptDir + '/TMSdata', 'r')
    TMSdata = TMSfile.readlines()
    TMSfile.close()

    for i, line in enumerate(TMSdata):
        buf = line.split(' ')
        if len(buf) > 1:
            if settings.Solvent != '':
                if buf[0].lower() == settings.nFunctional.lower() and \
                        buf[1].lower() == settings.nBasisSet.lower() and \
                        buf[2].lower() == settings.Solvent.lower():
                    print("Setting TMS references to " + buf[3] + " and " + \
                          buf[4] + "\n")
                    TMS_SC_C13 = float(buf[3])
                    TMS_SC_H1 = float(buf[4])
                    return TMS_SC_C13, TMS_SC_H1
            else:
                if buf[0].lower() == settings.nFunctional.lower() and \
                        buf[1].lower() == settings.nBasisSet.lower() and \
                        buf[2].lower() == 'none':
                    print("Setting TMS references to " + buf[3] + " and " + \
                          buf[4] + "\n")
                    TMS_SC_C13 = float(buf[3])
                    TMS_SC_H1 = float(buf[4])
                    return TMS_SC_C13, TMS_SC_H1

    print("No TMS reference data found for these conditions, using defaults\n")
    print("Unscaled shifts might be inaccurate, use of unscaled models is" + \
          " not recommended.")

    return settings.TMS_SC_C13, settings.TMS_SC_H1


def NMRDataValid(Isomers):

    for isomer in Isomers:
        if (len(isomer.ConformerShieldings) == 0):
            return False

    return True


def CalcNMRShifts(Isomers, settings):

    print('WARNING: NMR shift calculation currently ignores the instruction to exclude atoms from analysis')
    for i, iso in enumerate(Isomers):

        BShieldings = iso.BoltzmannShieldings
        Cvalues = []
        Hvalues = []

        for a, atom in enumerate(iso.Atoms):

            if atom == 'C':
                shift = (settings.TMS_SC_C13-BShieldings[a]) / (1-(settings.TMS_SC_C13/10**6))
                Cvalues.append(shift)

            if atom == 'H':
                shift = (settings.TMS_SC_H1-BShieldings[a]) / (1-(settings.TMS_SC_H1/10**6))
                Hvalues.append(shift)

        Isomers[i].Cshifts = Cvalues
        Isomers[i].Hshifts = Hvalues

        print('C shifts for isomer ' + str(i) + ": ")
        print(', '.join(['{0:.3f}'.format(x) for x in Isomers[i].Cshifts]))

        print('H shifts for isomer ' + str(i) + ": ")
        print(', '.join(['{0:.3f}'.format(x) for x in Isomers[i].Hshifts]))

    return Isomers


def PrintConformationData(AllSigConfs):
    # Make a list of populations and corresponding files for reporting
    # significant conformations
    """from operator import itemgetter
    ConfsPops = [list(x) for x in zip(args, populations)]
    ConfsPops.sort(key=itemgetter(1), reverse=True)
    totpop = 0
    i = 0
    while totpop < 0.8:
        totpop += ConfsPops[i][1]
        i += 1
    SigConfs = ConfsPops[:i]"""
    for Es, pops in zip(RelEs, populations):
        print('\nConformer relative energies (kJ/mol): ' + \
            ', '.join(["{:5.2f}".format(float(x)) for x in Es]))

        print('\nPopulations (%): ' + \
            ', '.join(["{:4.1f}".format(float(x)*100) for x in pops]))

    for i, SigConfs in enumerate(AllSigConfs):
        print("\nNumber of significant conformers for isomer "\
            + str(i+1) + ": " + str(len(SigConfs)) + "\n(pop, filename)")
        for conf in SigConfs:
            print("   " + format(conf[1]*100, "4.2f") + "%   " + conf[0])
        print('----------------')
        print("   " + format(100*sum([x[1] for x in SigConfs]), "4.2f") +\
            "%   in total")


def RemoveEquivalents(Noutp, equivs, OldCval, OldHval, OldClabels, OldHlabels):
    Cvalues = list(OldCval)
    Hvalues = list(OldHval)
    Clabels = list(OldClabels)
    Hlabels = list(OldHlabels)
    
    for eqAtoms in equivs:

        eqSums = [0.0]*Noutp
        eqAvgs = [0.0]*Noutp

        if eqAtoms[0][0] == 'H':
            #print eqAtoms, Hlabels
            for atom in eqAtoms:
                eqIndex = Hlabels.index(atom)
                for ds in range(0, Noutp):
                    eqSums[ds] = eqSums[ds] + Hvalues[ds][eqIndex]
            for ds in range(0, Noutp):
                eqAvgs[ds] = eqSums[ds]/len(eqAtoms)

            #Place the new average value in the first atom shifts place
            target_index = Hlabels.index(eqAtoms[0])
            for ds in range(0, Noutp):
                Hvalues[ds][target_index] = eqAvgs[ds]

            #Delete the redundant atoms from the computed list
            #start with second atom - e.g. don't delete the original one
            for atom in range(1, len(eqAtoms)):
                del_index = Hlabels.index(eqAtoms[atom])
                del Hlabels[del_index]
                for ds in range(0, Noutp):
                    del Hvalues[ds][del_index]

        if eqAtoms[0][0] == 'C':
            for atom in eqAtoms:
                eqIndex = Clabels.index(atom)
                for ds in range(0, Noutp):
                    eqSums[ds] = eqSums[ds] + Cvalues[ds][eqIndex]
            for ds in range(0, Noutp):
                eqAvgs[ds] = eqSums[ds]/len(eqAtoms)

            #Place the new average value in the first atom shifts place
            target_index = Clabels.index(eqAtoms[0])
            for ds in range(0, Noutp):
                Cvalues[ds][target_index] = eqAvgs[ds]

            #Delete the redundant atoms from the computed list
            #start with second atom - e.g. don't delete the original one
            for atom in range(1, len(eqAtoms)):
                del_index = Clabels.index(eqAtoms[atom])
                del Clabels[del_index]
                for ds in range(0, Noutp):
                    del Cvalues[ds][del_index]
                    
    return Cvalues, Hvalues, Clabels, Hlabels
    

def MAE(L1, L2):

    if len(L1) != len(L2):
        return -1
    else:
        L = []
        for i in range(0, len(L1)):
            L.append(abs(L1[i]-L2[i]))
        return sum(L)/len(L)


def RMSE(L1, L2):

    if len(L1) != len(L2):
        return -1
    else:
        L = []
        for i in range(0, len(L1)):
            L.append((L1[i]-L2[i])**2)
        return math.sqrt(sum(L)/len(L))
