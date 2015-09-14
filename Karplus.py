# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 15:20:18 2014

@author: ke291

Gets called by PyDP4.py if J value analysis based on Karplus equation is used
"""

import numpy as np
from math import sqrt, pi, cos, sin, acos
import sys
sys.path.append(
    '/home/ke291/Tools/openbabel-install/lib/python2.7/site-packages/')
from openbabel import *


def Karplus(f):

    """
    Scan molecule for sp3 protons (protons->carbons->sp3)
    Check if they have dihedral sp3 protons
    calculate dihedral matrix
    calculate J value matrix based on karplus equation

    """
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()

    obconversion.ReadFile(obmol, f)

    obmol.ConnectTheDots()
    obmol.Kekulize()
    
    DihedralHs = []
    for atom in OBMolAtomIter(obmol):
        if atom.GetAtomicNum() == 1:
            DihedNeighbours = GetDihedralHs(atom)
            if DihedNeighbours!=0:
                DihedralHs.append([atom.GetIdx()] + DihedNeighbours)                
    
    print DihedralHs
    
    DihedralMatrix, labels = CalcDihedralMatrix(obmol, DihedralHs)
    
    print labels
    print DihedralMatrix
    """
    Jmatrix, Jlabels = CalculateJs(DihedralMatrix)
    
    return Jmatrix, Jlabels
"""

def GetDihedralHs(atom):
    sp3carbonFound = False
    #find the proton neigbour, check that it is sp3 carbon
    for NbrAtom in OBAtomAtomIter(atom):
        if NbrAtom.GetAtomicNum() == 6 and NbrAtom.GetHyb() == 3:
            carbon = NbrAtom
            sp3carbonFound = True
            break
    
    if not sp3carbonFound:
        return 0
    #find the carbon neibours, check that they are sp3 hybridized
    ncarbons = []
    for NbrAtom in OBAtomAtomIter(carbon):
        if NbrAtom.GetAtomicNum() == 6 and NbrAtom.GetHyb() == 3:
            ncarbons.append(NbrAtom)
    #check that the carbon neighbours have protons and return them
    nHs = []
    for ncarbon in ncarbons:
        for NbrAtom in OBAtomAtomIter(ncarbon):
            if NbrAtom.GetAtomicNum() == 1:
                nHs.append(NbrAtom.GetIdx())
    
    if len(nHs) == 0:
        return 0
    else:
        return nHs


def CalcDihedralMatrix(obmol, DihedralHs):
    mat = [[0.0 for x in range(len(DihedralHs))] for y in range(len(DihedralHs))]
    labels = [x[0] for x in DihedralHs]
    
    for x, row in enumerate(DihedralHs):
        atom1 = row[0]
        for atom2 in row[1:]:
            y = labels.index(atom2)
            #angle = CalcDihedral(atom1, atom2)
            angle = 57.2957795*CalcDihedral(obmol.GetAtom(atom1),
                                            obmol.GetAtom(atom2))
            mat[x][y] = angle
            mat[y][x] = angle
    
    return mat, labels


def CalcDihedral(atom1, atom2):
    
    #Get the 2 carbon atoms joining the protons
    for NbrAtom in OBAtomAtomIter(atom1):
        if NbrAtom.GetAtomicNum() == 6 and NbrAtom.GetHyb() == 3:
            carbon1 = NbrAtom
            break
    for NbrAtom in OBAtomAtomIter(atom2):
        if NbrAtom.GetAtomicNum() == 6 and NbrAtom.GetHyb() == 3:
            carbon2 = NbrAtom
            break
    #Get the 2 planes formed by the carbons and each of the protons
    normal1, d = FindPlane(atom1, carbon1, carbon2)
    normal2, d = FindPlane(atom2, carbon1, carbon2)
    
    return AngleBetweenPlanes(normal1, normal2)


def AngleBetweenPlanes(normal1, normal2):
    
    dotprod = sum((a*b) for a, b in zip(normal1, normal2))
    length1 = sqrt(sum([x**2 for x in normal1]))
    length2 = sqrt(sum([x**2 for x in normal2]))
    
    return acos(dotprod/(length1*length2))


def GetUnitVector(Atom1, Atom2):
    vector = []
    vector.append(Atom2.x() - Atom1.x())
    vector.append(Atom2.y() - Atom1.y())
    vector.append(Atom2.z() - Atom1.z())

    length = np.linalg.norm(vector)
    return [x/length for x in vector]


#Given 3 atoms, finds a plane defined by a normal vector and d
def FindPlane(atom1, atom2, atom3):

    vector1 = [atom2.x() - atom1.x(), atom2.y() - atom1.y(),
               atom2.z() - atom1.z()]
    vector2 = [atom3.x() - atom1.x(), atom3.y() - atom1.y(),
               atom3.z() - atom1.z()]
    cross_product = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
                     -1 * vector1[0] * vector2[2] - vector1[2] * vector2[0],
                     vector1[0] * vector2[1] - vector1[1] * vector2[0]]

    d = cross_product[0] * atom1.x() - cross_product[1] * atom1.y() + \
        cross_product[2] * atom1.z()

    return cross_product, d


def crossproduct(v1, v2):
    product = [0, 0, 0]
    product[0] = v1[1]*v2[2]-v1[2]*v2[1]
    product[1] = v1[2]*v2[0]-v1[0]*v2[2]
    product[2] = v1[0]*v2[1]-v1[1]*v2[0]
    return product


def dotproduct(v1, v2):
        return sum((a*b) for a, b in zip(v1, v2))


def length(v):
    return sqrt(dotproduct(v, v))


def angle(v1, v2):
    return acos(dotproduct(v1, v2) / (length(v1) * length(v2)))