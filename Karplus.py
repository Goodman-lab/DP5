# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 15:20:18 2014

@author: ke291

Gets called by PyDP4.py if J value analysis based on Karplus equation is used
Based on C. A. G. Haasnoot, F. A. A. M. de Leeuw and C. Altona;
Tetrahedron 36 2783-2792 (1980)
"""

import numpy as np
from math import sqrt, pi, cos, sin, acos, atan2
import sys

try:
    from openbabel.openbabel import OBConversion, OBMol, OBAtomAtomIter, OBMolAtomIter
except ImportError:
    from openbabel import *

#Parameter sets B-E for generalised Karplus equation from
#C. A. G. Haasnoot, F. A. A. M. de Leeuw and C. Altona
#Tetrahedron 36 2783-2792 (1980)
Params = [[13.70,-0.73,0.0,0.56,-2.47,16.9,0.14],
          [13.89,-0.98,0.0,1.02,-3.40,14.9,0.24],
          [13.22,-0.99,0.0,0.87,-2.46,19.9,0.00],
          [13.24,-0.91,0.0,0.53,-2.41,15.5,0.19]]

#Electronegativities up to iodine
ENs = [2.20,0.00,0.98,1.57,2.04,2.55,3.04,3.44,3.98,0.00,
       0.93,1.31,1.61,1.90,2.19,2.58,3.16,0.00,0.82,1.00,
       1.36,1.54,1.63,1.66,1.55,1.83,1.88,1.91,1.90,1.65,
       1.81,2.01,2.18,2.55,2.96,3.00,0.82,0.95,1.22,1.33,
       1.60,2.16,1.90,2.20,2.28,2.20,1.93,1.69,1.78,1.96,
       2.05,2.10,2.66]

#input format - "g09" for gaussian, "nwo" for NWChem
def Karplus(f, inputformat):

    obconversion = OBConversion()
    obconversion.SetInFormat(inputformat)
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
    
    if len(DihedralHs)==0:
        print("No dihedral protons found, Karplus J value prediction " + \
            "impossible, quitting.")
        quit()
    
    Jmatrix, Jlabels = CalcJMatrix(obmol, DihedralHs)
    Jlabels = [str(x) for x in Jlabels]
    return Jmatrix, Jlabels


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


def CalcJMatrix(obmol, DihedralHs):
    mat = [[0.0 for x in range(len(DihedralHs))] for y in range(len(DihedralHs))]
    labels = [x[0] for x in DihedralHs]
    
    for x, row in enumerate(DihedralHs):
        atom1 = row[0]
        for atom2 in row[1:]:
            y = labels.index(atom2)
            J = CalcJ(obmol.GetAtom(atom1), obmol.GetAtom(atom2))
            mat[x][y] = J
            mat[y][x] = J
    
    return mat, labels


def CalcJ(atom1, atom2):
    
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
    normal1, d = FindPlane(carbon1, atom1, carbon2)
    normal2, d = FindPlane(carbon2, atom2, carbon1)
    
    sign = VectAngleSign(BondVect(carbon1, atom1), BondVect(carbon2, atom2),
                         BondVect(carbon1, carbon2))
    angle = VectorAngle(normal1, normal2)
    if dotproduct(BondVect(carbon1, atom1), BondVect(carbon2, atom2)) > 0 and \
        angle < 0.5*pi:
        dihedral = angle*sign
    else:
        dihedral = (pi-angle)*sign
    
    #Count the substituents, get their electronegativities
    nSubst = 0
    relENs = []
    ENcorrs = []
    signs = []
    
    for NbrAtom in OBAtomAtomIter(carbon1):
        ANum = NbrAtom.GetAtomicNum()
        if NbrAtom.GetIdx() != carbon2.GetIdx() and ANum != 1:
            nSubst += 1
            relEN, ENcorr = GetENCorrection(NbrAtom, carbon1)
            relENs.append(relEN)
            ENcorrs.append(ENcorr)
            signs.append(VectAngleSign(BondVect(carbon1, atom1),
                                       BondVect(carbon1, NbrAtom),
                                       BondVect(carbon1, carbon2)))
                        
    for NbrAtom in OBAtomAtomIter(carbon2):
        ANum = NbrAtom.GetAtomicNum()
        if NbrAtom.GetIdx() != carbon1.GetIdx() and ANum != 1:
            nSubst += 1
            relEN, ENcorr = GetENCorrection(NbrAtom, carbon2)
            relENs.append(relEN)
            ENcorrs.append(ENcorr)
            signs.append(VectAngleSign(BondVect(carbon2, atom2),
                                       BondVect(carbon2, NbrAtom),
                                       BondVect(carbon2, carbon1)))
    
    SubstEffects = 0
    if nSubst < 2:
        for relEN, ENcorr, sign in zip(relENs, ENcorrs, signs):
            corrRelEN = relEN - (Params[0][6])*ENcorr
            SubstEffects += SubstEffect(dihedral, corrRelEN, sign, 0)
            
        J = Params[0][0]*cos(dihedral)**2 + Params[0][1]*cos(dihedral)
    else:
        for relEN, sign in zip(relENs, signs):
            corrRelEN = relEN - (Params[nSubst-1][6])*ENcorr
            SubstEffects += SubstEffect(dihedral, corrRelEN, sign, nSubst-1)
            
        J = Params[nSubst-1][0]*cos(dihedral)**2 + Params[nSubst-1][1]*cos(dihedral)
    
    return J + SubstEffects


def GetENCorrection(atom, satom):
    Anum = atom.GetAtomicNum()
    if Anum != 6:
        return ENs[Anum-1]-ENs[0], 0.0
    else:
        ENcorr = 0.0
        for NbrAtom in OBAtomAtomIter(atom):
            if NbrAtom.GetIdx() != satom.GetIdx():
                neighbAnum = NbrAtom.GetAtomicNum()
                ENcorr += ENs[neighbAnum-1]-ENs[0]
        
        return ENs[Anum-1]-ENs[0], ENcorr


def SubstEffect(dihedral, relEN, sign, pset):
    
    return relEN*(Params[pset][3] + Params[pset][4]*
                  cos(sign*dihedral + Deg2Rad(Params[pset][5])*abs(relEN))**2)


def BondVect(atom1, atom2):
    return [atom2.x() - atom1.x(),atom2.y() - atom1.y(),atom2.z() - atom1.z()]


def VectorAngle(v1, v2):
    
    cp = crossproduct(v1, v2)
    return atan2(length(cp),dotproduct(v1, v2))


def VectorAngle2(normal1, normal2):
    
    dotprod = sum((a*b) for a, b in zip(normal1, normal2))
    length1 = sqrt(sum([x**2 for x in normal1]))
    length2 = sqrt(sum([x**2 for x in normal2]))
    
    return acos(dotprod/(length1*length2))


def VectAngleSign(v1, v2, axis):
    
    crossprod = crossproduct(v1, v2)
    return np.sign(dotproduct(crossprod, axis))


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
    cross_product = crossproduct(vector1, vector2)
    d = cross_product[0] * atom1.x() - cross_product[1] * atom1.y() + \
        cross_product[2] * atom1.z()
    
    return cross_product, d


def Deg2Rad(x):
    return x*0.0174532925


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