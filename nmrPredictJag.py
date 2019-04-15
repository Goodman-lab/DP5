# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 12:46:34 2014

@author: ke291

Extracts NMR shifts from Gaussian output files
"""

import math

gasConstant = 8.3145
temperature = 298.15
hartreeEnergy = 2625.499629554010

def main(settings, *args):

    labels = []
    allshieldings = []
    energies = []

    #For every file run ReadShieldings and store the extracted shielding
    #constants
    for f in args:
        (labels, shieldings) = ReadShieldings(f)
        allshieldings.append(shieldings)
        energy = ReadEnergy(settings.EnergyDir + f)
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
    
    #Make a list of populations and corresponding files for reporting
    #significant conformations
    from operator import itemgetter
    ConfsPops = [list(x) for x in zip(args, populations)]
    ConfsPops.sort(key=itemgetter(1), reverse=True)
    totpop = 0
    i = 0
    while totpop<0.8:
        totpop += ConfsPops[i][1]
        i += 1
    SigConfs = ConfsPops[:i]
        
    #Calculate Boltzmann weighed shielding constants
    #by summing the shifts multiplied by the isomers population
    Natoms = len(labels)
    BoltzmannShieldings = []

    for atom in range(Natoms):
        shielding = 0
        for conformer in range(len(allshieldings)):
            shielding = shielding + allshieldings[conformer][atom] * \
                populations[conformer]
        BoltzmannShieldings.append(shielding)

    return (relEs, populations, labels, BoltzmannShieldings, SigConfs)


def ReadEnergy(JOutpFile):
    jagfile = open(JOutpFile + '.out', 'r')
    JOutp = jagfile.readlines()
    jagfile.close()
    energies = []
    
    for i, line in enumerate(JOutp):
        if 'SCFE:  Solution phase energy: DFT' in line:
            start = line.index(')  ')
            end = line.index(' hartrees')
            energies.append(float(line[start+4:end]))
    
    return energies[-1]


def ReadShieldings(JOutpFile):
    
    print(JOutpFile)
    jagfile = open(JOutpFile + '.out', 'r')
    JOutp = jagfile.readlines()

    index = 0
    shieldings = []
    labels = []
    #Find the NMR shielding calculation section
    while not 'NMR Properties' in JOutp[index]:
        index = index + 1
    #Read shielding constants and labels
    for i in range(index, len(JOutp)):
        line = JOutp[i]
        if 'NMR Properties' in line:
            data = [_f for _f in line.split(' ') if _f]
            labels.append(data[-2])
            line = JOutp[i+3]
            data = [_f for _f in line.split(' ') if _f]
            shieldings.append(float(data[-1]))
            
    return labels, shieldings
