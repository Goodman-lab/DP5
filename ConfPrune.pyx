# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 12:33:18 2015

@author: ke291

Cython file for conformer alignment and RMSD pruning. Called by Gaussian.py
and NWChem.py
"""

from math import sqrt
from libc.stdlib cimport malloc, free
"""
Main function of this file. Takes an array of [x, y, z] for each molecule,
as well as a list of atom symbols in the same order as their coordinates

Translates both molecules to origin, gets the rotation matrix from quatfit
and rotates the second molecule to align with the first
"""
def RMSDPrune(conformers, atoms, cutoff):
    
    ToDel = []
    cdef long int c1, c2, c
    cdef long int l = len(conformers)
    cdef double res
    cdef double *RMSDMatrix = <double *>malloc(l * l * sizeof(double))
    
    for c1 in range(0, l):
        for c2 in range(c1, l):
            if c1==c2:
                RMSDMatrix[c2 + c1*l]=0.0
            else:
                res = AlignedRMS(conformers[c1], conformers[c2], atoms)
                RMSDMatrix[c2 + c1*l] = res
                RMSDMatrix[c1 + c2*l] = res
    #Check for similar conformations
    for c1 in range(0, l):
        for c2 in range(0, l):
            if c1!=c2 and (not c1 in ToDel) and (not c2 in ToDel):
                #if Align.AlignedRMS(conformers[c1], conformers[c2], atoms) < cutoff:
                #    ToDel.append(c2)
                if RMSDMatrix[c2 + c1*l]<cutoff:
                    ToDel.append(c2)
    
    #Compose set of non-redundant conformations
    PrunedConformers = []
    for c in range(0, l):
        if not c in ToDel:
            PrunedConformers.append(conformers[c])
            
    #return PrunedConformers
    free(RMSDMatrix)
    return PrunedConformers

def AdaptRMSDPrune(conformers, atoms, cutoff, ConfLimit):
    
    ToDel = []
    cdef long int c1, c2
    cdef long int l = len(conformers)
    cdef double res
    cdef double *RMSDMatrix = <double *>malloc(l * l * sizeof(double))
    
    for c1 in range(0, l):
        for c2 in range(c1, l):
            if c1==c2:
                RMSDMatrix[c2 + c1*l]=0.0
            else:
                res = AlignedRMS(conformers[c1], conformers[c2], atoms)
                RMSDMatrix[c2 + c1*l] = res
                RMSDMatrix[c1 + c2*l] = res
    #Check for similar conformations
    for c1 in range(0, l):
        for c2 in range(0, l):
            if c1!=c2 and (not c1 in ToDel) and (not c2 in ToDel):
                if RMSDMatrix[c2 + c1*l]<cutoff:
                    ToDel.append(c2)
    AdjCutoff = cutoff
    while (l-len(ToDel))>ConfLimit:
        AdjCutoff +=0.2
        ToDel = []
        for c1 in range(0, l):
            for c2 in range(0, l):
                if c1!=c2 and (not c1 in ToDel) and (not c2 in ToDel):
                    if RMSDMatrix[c2 + c1*l]<AdjCutoff:
                        ToDel.append(c2)
              
    #return PrunedConformers
    free(RMSDMatrix)
    return AdjCutoff


def StrictRMSDPrune(conformers, atoms, cutoff, ConfLimit):
    
    ToDel = []
    cdef long int c1, c2
    cdef long int l = len(conformers)
    cdef double res
    cdef double *RMSDMatrix = <double *>malloc(l * l * sizeof(double))
    
    for c1 in range(0, l):
        for c2 in range(c1, l):
            if c1==c2:
                RMSDMatrix[c2 + c1*l]=0.0
            else:
                res = AlignedRMS(conformers[c1], conformers[c2], atoms)
                RMSDMatrix[c2 + c1*l] = res
                RMSDMatrix[c1 + c2*l] = res
    #Check for similar conformations
    for c1 in range(0, l):
        for c2 in range(0, l):
            if c1!=c2 and (not c1 in ToDel) and (not c2 in ToDel):
                if RMSDMatrix[c2 + c1*l]<cutoff:
                    ToDel.append(c2)
    AdjCutoff = cutoff
    while (l-len(ToDel))>ConfLimit:
        AdjCutoff +=0.2
        ToDel = []
        for c1 in range(0, l):
            for c2 in range(0, l):
                if c1!=c2 and (not c1 in ToDel) and (not c2 in ToDel):
                    if RMSDMatrix[c2 + c1*l]<AdjCutoff:
                        ToDel.append(c2)
    #Compose set of non-redundant conformations
    PrunedConformers = []
    for c in range(0, l):
        if not c in ToDel:
            PrunedConformers.append(conformers[c])
            
    #return PrunedConformers
    free(RMSDMatrix)
    return PrunedConformers, AdjCutoff


def AlignMolecules(mol1, mol2, atoms):
    
    w= []
    #prepare weights
    for a in atoms:
        atomnum = GetAtomNum(a)
        #w.append(GetAtomWeight(atomnum))
        w.append(1)
    
    #move both molecules to origin
    mol1 = Move2Origin(mol1, w)
    mol2 = Move2Origin(mol2, w)
    
    #QuatrFit to get the rotation matrix
    u = qtrfit(mol1, mol2, w)
    
    #RotateMolecule
    mol1 = RotMol(mol1, u)
    
    return mol1, mol2, w

#Converts to pure geometry data, aligns molecules and
#returns the RMSD of the aligned molecules
def AlignedRMS(mol1, mol2, atoms):
    
    #Convert to pure geometry data from [atomsym,text x,text y,text z]
    mol1b = []
    mol2b = []
    for a in mol1:
        mol1b.append([float(x) for x in a[1:]])
        
    for a in mol2:
        mol2b.append([float(x) for x in a[1:]])
    
    (amol1, amol2, w) = AlignMolecules(mol1b, mol2b, atoms)
    
    return RMSMol(amol1, amol2, w)

#calculates the RMS between the two molecules
def RMSMol(mol1, mol2, w):
    e = [0.0, 0.0, 0.0]
    
    for i in range(0, len(mol1)):
        e[0] += w[i] * (mol1[i][0] - mol2[i][0])**2
        e[1] += w[i] * (mol1[i][1] - mol2[i][1])**2
        e[2] += w[i] * (mol1[i][2] - mol2[i][2])**2
    
    return sqrt(sum(e)/sum(w))

def NMR_RMS(s1, s2):
    
    e = 0
    
    for i in range(0, len(s1)):
        e += (s1[i] - s2[i])**2

    return sqrt(e/len(s1))
#Rotates molecule, given an array of atoms and a rotation matrix
#molecule is an array of form a[atom][coordinate]

def RotMol(mol, u):
    
    t = [0.0, 0.0, 0.0]
    
    for i in range(0, len(mol)):
        t[0] = u[0][0] * mol[i][0] + u[0][1] * mol[i][1] + u[0][2] * mol[i][2]
        t[1] = u[1][0] * mol[i][0] + u[1][1] * mol[i][1] + u[1][2] * mol[i][2]
        t[2] = u[2][0] * mol[i][0] + u[2][1] * mol[i][1] + u[2][2] * mol[i][2]
        mol[i][0] = t[0]
        mol[i][1] = t[1]
        mol[i][2] = t[2]
    
    return mol

#Translates the molecule so that it's centroid is at the origin
#weights can be atomic weights or numbers?
def Move2Origin(mol, weights):
    
    xsum = 0.0
    ysum = 0.0
    zsum = 0.0
    wsum = 0.0
    
    for i in range(0, len(mol)):
        xsum += mol[i][0] * sqrt(weights[i])
        ysum += mol[i][1] * sqrt(weights[i])
        zsum += mol[i][2] * sqrt(weights[i])
        wsum += sqrt(weights[i])
    
    xsum /= wsum
    ysum /= wsum
    zsum /= wsum
    
    for i in range(0, len(mol)):
        mol[i][0] -= xsum
        mol[i][1] -= ysum
        mol[i][2] -= zsum
    
    return mol

def GetAtomNum(symbol):
    
    Lookup = ['H', 'He', 'Li', 'Be', 'B', 'C','N','O','F','Ne','Na','Mg','Al',\
            'Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co',\
            'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr',\
            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I',\
            'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',\
            'Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au',\
            'Hg','Tl','Pb','Bi','Po','At','Rn']
    
    if symbol in Lookup:
        return Lookup.index(symbol) + 1
    else:
        return 0

def GetAtomWeight(AtomNum):
    
    Lookup = [1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35,\
            40, 39, 40, 45,  48, 51, 52, 55, 56, 59, 59, 64, 65, 70, 73, 75,\
            79, 80, 84, 85, 88, 89, 91, 93, 96, 98, 101, 103, 106, 108, 112,\
            115, 119, 122, 128, 127, 131, 133, 137]
    return Lookup[AtomNum-1]

"""
 * This routine is a python translation of the QTRFIT Fortran routine
 * of David Heisterberg
 *
 * David J. Heisterberg
 * The Ohio Supercomputer Center, 1224 Kinnear Rd,
 * Columbus, OH  43212-1163
 * (614)292-6036; djh@osc.edu
 *
 * Python translation by K Ermanis 2015
 *
 * QTRFIT
 * Find the quaternion, q,[and left rotation matrix, u] 
 * that minimizes
 *   |qTBq - A| ^ 2   [|uB - A| ^ 2]
 * This is equivalent to maximizing Re (qTBTqA).
 * This is equivalent to finding the largest 
 * eigenvalue and corresponding
 * eigenvector of the matrix
 *  [V2   VUx  VUy  VUz ]
 *  [VUx  Ux2  UxUy UzUx]
 *  [VUy  UxUy Uy2  UyUz]
 *  [VUz  UzUx UyUz Uz2 ]
 * where
 *   V2   = Bx Ax + By Ay + Bz Az
 *   Ux2  = Bx Ax - By Ay - Bz Az
 *   Uy2  = By Ay - Bz Az - Bx Ax
 *   Uz2  = Bz Az - Bx Ax - By Ay
 *   VUx  = Bz Ay - By Az
 *   VUy  = Bx Az - Bz Ax
 *   VUz  = By Ax - Bx Ay
 *   UxUy = Bx Ay + By Ax
 *   UyUz = By Az + Bz Ay
 *   UzUx = Bz Ax + Bx Az
 * The left rotation matrix, u, is obtained from q by
 *   u = qT1q
 * INPUT
 *   n      - number of points
 *   b      - molecule to be rotated
 *   a      - reference molecule
 *   w      - weight vector
 * OUTPUT
 *   z      - eigenvalue
 *   q      - the best-fit quaternion
 *   u      - the best-fit left rotation matrix
 *   nr     - number of jacobi sweeps required
 
 KE:
 structure of molecule data:
 a[atomnr][coordinate]
 """
cdef double v[4][4]
cdef double d[4]
    
def qtrfit(mol, refmol, w):
    
    cdef int n = len(mol)
    cdef double xxyx = 0.0, xxyy = 0.0, xxyz = 0.0, xyyx = 0.0, xyyy = 0.0
    cdef double xyyz = 0.0, xzyx = 0.0, xzyy = 0.0, xzyz = 0.0
    cdef double c[4][4]
    cdef double q[4]
    cdef int i, j
    #generate the upper triangle of the quadratic form matrix:
    for i in range(0,n):
        xxyx += mol[i][0] * refmol[i][0] * w[i]
        xxyy += mol[i][0] * refmol[i][1] * w[i]
        xxyz += mol[i][0] * refmol[i][2] * w[i]
        xyyx += mol[i][1] * refmol[i][0] * w[i]
        xyyy += mol[i][1] * refmol[i][1] * w[i]
        xyyz += mol[i][1] * refmol[i][2] * w[i]
        xzyx += mol[i][2] * refmol[i][0] * w[i]
        xzyy += mol[i][2] * refmol[i][1] * w[i]
        xzyz += mol[i][2] * refmol[i][2] * w[i]
    
    c[0][0] = xxyx + xyyy + xzyz
    c[0][1] = xzyy - xyyz
    c[1][1] = xxyx - xyyy - xzyz
    c[0][2] = xxyz - xzyx
    c[1][2] = xxyy + xyyx
    c[2][2] = xyyy - xzyz - xxyx
    c[0][3] = xyyx - xxyy
    c[1][3] = xzyx + xxyz
    c[2][3] = xyyz + xzyy
    c[3][3] = xzyz - xxyx - xyyy
    
    #Diagonalize c
    cdef double **ret
    ret = qtrjac(c, 4, 4)
        
    #extract the desired quaternion
    cdef double z = d[3]
    q[0] = v[0][3]
    q[1] = v[1][3]
    q[2] = v[2][3]
    q[3] = v[3][3]
    
    #generate the rotation matrix
    u = q2mat(q)
    
    return u

cdef double** qtrjac(double a[4][4], int n, int np):
    
    """v = [[0.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0],\
        [0.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0]]
    d = [0.0, 0.0, 0.0, 0.0]
    onorm = 0.0
    dnorm = 0.0
    b = 0.0
    dma = 0.0
    q = 0.0
    t= 0.0
    c = 0.0
    s = 0.0
    atemp = 0.0
    vtemp = 0.0
    dtemp = 0.0
    nrot = 128"""
    cdef double onorm, dnorm, b, dma, q, t, c, s, atemp, vtemp, dtemp
    cdef int nrot = 128
    cdef int j, i, l, k
    for j in range(0,n):
        for i in range(0,n):
            v[i][j]=0.0
        v[j][j] = 1.0
        d[j] = a[j][j]

    for l in range(0, nrot):
        dnorm = 0.0
        onorm = 0.0
        for j in range(0, n):
            dnorm += abs(d[j])
            for i in range(j-1):
                onorm += abs(a[i][j])
        if (onorm + dnorm) <= dnorm:
            break
        
        for j in range(1, n):
            
            for i in range(0,j-1):
                
                b = a[i][j]
                if (abs(b)>0.0):
                    dma = d[j] - d[i]
                    if (abs(dma)+abs(b) < abs(dma)):
                        t = b/dma
                    else:
                        q = 0.5 * dma / b
                        t = fsign(1.0 / (abs(q) + sqrt(1.0 + q*q)), q)
                    
                    c = 1.0 / sqrt(t*t + 1.0)
                    s = t * c
                    a[i][j] = 0.0
                    
                    for k in range(0, i-1):
                        atemp = c * a[k][i] - s * a[k][j]
                        a[k][j] = s * a[k][i] + c * a[k][j]
                        a[k][i] = atemp
                        
                    for k in range(i+1, j):
                        atemp = c * a[i][k] - s * a[k][j]
                        a[k][j] = s * a[i][k] + c * a[k][j]
                        a[i][k] = atemp
                        
                    for k in range(j, n):
                        atemp = c * a[i][k] - s * a[j][k]
                        a[j][k] = s * a[i][k] + c * a[j][k]
                        a[i][k] = atemp
                    
                    for k in range(0, n):
                        vtemp = c * v[k][i] - s * v[k][j]
                        v[k][j] = s * v[k][i] + c * v[k][j]
                        v[k][i] = vtemp
                    
                    dtemp = c * c * d[i] + s * s * d[j] - 2.0 + c * s * b
                    d[j] = s * s * d[i] + c * c * d[j] + 2.0 * c * s * b
                    d[i] = dtemp
    nrot = l
    
    for j in range(0, n-1):
        k = j
        dtemp = d[k]
        for i in range(j, n):
            if d[i]<dtemp:
                k = i
                dtemp = d[k]
        if k>j:
            d[k] = d[j]
            d[j] = dtemp
            for i in range(0, n):
                dtemp = v[i][k]
                v[i][k] = v[i][j]
                v[i][j] = dtemp
    
    cdef double ret[5][4]
    for i in range(0,4):
        for j in range(0,4):
            ret[i][j]=v[i][j]
    for i in range(0,4):
        ret[4][i] = d[i]
    

cdef q2mat(double q[4]):
    
    u = [[0.0, 0.0, 0.0],[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]]
    u[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3]
    u[1][0] = 2.0 *(q[1] * q[2] - q[0] * q[3])
    u[2][0] = 2.0 *(q[1] * q[3] + q[0] * q[2])
    u[0][1] = 2.0 *(q[2] * q[1] + q[0] * q[3])
    u[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3]
    u[2][1] = 2.0 *(q[2] * q[3] - q[0] * q[1])
    u[0][2] = 2.0 *(q[3] * q[1] - q[0] * q[2])
    u[1][2] = 2.0 *(q[3] * q[2] + q[0] * q[1])
    u[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]
    
    return u

def fsign(a, b):
    if (b>=0):
        return a
    else:
        return -a
