# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 12:46:34 2014

@author: ke291

Extracts NMR shifts from NWChem output files
"""

import sys
import math

gasConstant = 8.3145
temperature = 298.15
hartreeEnergy = 2625.499629554010

def main(*args):

    labels = []
    allshieldings = []
    energies = []

    #For every file run ReadShieldings and store the extracted shielding
    #constants
    for f in args:
        (labels, shieldings, energy) = ReadShieldings(f)
        allshieldings.append(shieldings)
        energies.append(energy)

    minE = min(energies)

    relEs = []

    #Calculate rel. energies in kJ/mol
    for e in energies:
        relEs.append(((e-minE)*hartreeEnergy))

    populations = []

    #Calculate Boltzmann populations
    for e in relEs:
        populations.append(math.exp(-e*1000/(gasConstant*temperature)))

    q = sum(populations)

    for i in range(0, len(populations)):
        populations[i] = populations[i]/q

    #Calculate Boltzmann weighed shielding constants
    #by summing the shifts multiplied by the isomers population
    Natoms = len(labels)
    BoltzmannShieldings = []

    for atom in range(0, Natoms):
        shielding = 0
        for conformer in range(0, len(allshieldings)):
            shielding = shielding + allshieldings[conformer][atom] * \
                populations[conformer]
        BoltzmannShieldings.append(shielding)

    return (relEs, populations, labels, BoltzmannShieldings)


def ReadShieldings(NWOutpFile):

    nwfile = open(NWOutpFile + '.nwo', 'r')
    NWOutp = nwfile.readlines()

    index = 0
    shieldings = []
    labels = []
    energies = []

    #Find the NMR shielding calculation section
    while not 'Chemical Shielding' in NWOutp[index]:
        index = index + 1

    for line in NWOutp:
        if 'Total DFT energy' in line:
            start = line.index('Total')
            energies.append(float(line[start+19:]))

    #Read shielding constants and labels
    for line in NWOutp[index:]:
        if 'isotropic' in line:
            start = line.index('isotropic')
            shieldings.append(float(line[start+13:]))
        if 'Atom' in line:
            start = line.index('Atom')
            labels.append(line[start+12] + line[start+7:start+10].strip())

    nwfile.close()

    return labels, shieldings, energies[-1]

if __name__ == '__main__':
    #print sys.argv
    cpArgs = sys.argv[1:]
    main(*cpArgs)
