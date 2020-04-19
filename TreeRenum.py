#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:54:17 2015

@author: ke291
"""
import sys

def FindAllPaths(molgraph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return [path]
    if start not in [x[0] for x in molgraph]:
        print("No such node in graph")
        return []
    paths = []
    for node in molgraph[start-1][1:]:
        if node not in path:
            newpaths = FindAllPaths(molgraph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths


def FindTerminatingPaths(molgraph, start, trunk, path=[]):
    path = path + [start]
    if start not in [x[0] for x in molgraph]:
        print("No such node in graph")
        return []
    paths = []
    for node in molgraph[start-1][1:]:
        if node not in path and node not in trunk:
            newpaths = FindTerminatingPaths(molgraph, node, trunk, path)
            for newpath in newpaths:
                paths.append(newpath)
    if paths == []:
        return [path]
    else:
        return paths


def FindTreeMap(f):
    #get molecular graph
    molgraph = GenMolGraph(f)

    #Find potential endpoints - atoms with only one neighbour
    endpoints = []
    for node in molgraph:
        if len(node) < 3:
            endpoints.append(node[0])

    #get the longest paths for all endpoint combinations
    maxpaths = []
    for i in range(0, len(endpoints)):
        for j in range(0, len(endpoints)):
            if i != j:
                maxpaths.append(max(FindAllPaths(molgraph, endpoints[i],
                                                 endpoints[j]), key=len))
    #get the absolute longest path possible in the molecule
    molmap = max(maxpaths, key=len)
    #Add longest branches to the longest trunk
    for atom in molmap:
        for node in molgraph[atom-1]:
            if node not in molmap:
                maxbranch = []
                branches = FindTerminatingPaths(molgraph, node, molmap)
                if branches != []:
                    maxbranch = max(branches, key=len)
                if maxbranch != []:
                    molmap.extend(maxbranch)
    return molmap


def TreeRenumSDF(f, ExpNMR):

    molmap = FindTreeMap(f)

    #Read molecule from file
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, f)

    temp = []
    anums = []
    for atom in OBMolAtomIter(obmol):
        temp.append(atom)
        anums.append(atom.GetAtomicNum())

    newmol = OBMol()
    for atom in molmap:
        newmol.AddAtom(temp[atom-1])

    #Restore the bonds
    newmol.ConnectTheDots()
    newmol.PerceiveBondOrders()
    #Write renumbered molecule to file
    obconversion.SetOutFormat("sdf")
    obconversion.WriteFile(newmol, f[:-4] + 'r.sdf')

    #Prepare NMR translation
    NMRmap = []
    i = 1
    for atom in molmap:
        anum = anums[atom-1]
        if anum == 1:
            NMRmap.append(['H' + str(atom), 'H' + str(i)])
        if anum == 6:
            NMRmap.append(['C' + str(atom), 'C' + str(i)])
        i += 1
    print(NMRmap)
    RenumNMR(ExpNMR, NMRmap)


def RenumNMR(ExpNMR, NMRmap):
    f = open(ExpNMR, 'r')
    NMRfile = f.read(1000000)
    f.close()

    print('\nOld NMR file:\n' + NMRfile)

    #Replace old atom labels with new atom labels
    #tag replacements with '_' to avoid double replacement
    for atom in NMRmap:
        NMRfile = NMRfile.replace(atom[0] + ')', atom[1] + '_)')
        NMRfile = NMRfile.replace(atom[0] + ' ', atom[1] + '_ ')
        NMRfile = NMRfile.replace(atom[0] + ',', atom[1] + '_,')
        NMRfile = NMRfile.replace(atom[0] + '\n', atom[1] + '_\n')

    #Strip temporary udnerscore tags
    NMRfile = NMRfile.replace('_', '')

    print('\nNew NMR file:\n' + NMRfile)
    f = open(ExpNMR + 'r', 'w')
    f.write(NMRfile)
    f.close()


def GenMolGraph(f):
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, f)

    molgraph = []

    for atom in OBMolAtomIter(obmol):
        idx = atom.GetIdx()
        molgraph.append([])
        molgraph[idx-1].append(idx)

        for NbrAtom in OBAtomAtomIter(atom):
            molgraph[idx-1].append(NbrAtom.GetIdx())

    return molgraph

if __name__ == '__main__':
    #print sys.argv
    ExpNMR = sys.argv[2]
    filename = sys.argv[1] + '.sdf'
    TreeRenumSDF(filename, ExpNMR)
    #main()