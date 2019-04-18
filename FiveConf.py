# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 15:20:18 2014

@author: ke291

Gets called by PyDP4.py if automatic 5-membered cycle corner-flipping is used.
"""

import numpy as np
from math import sqrt, pi, cos, sin, acos
import scipy.optimize as sciopt

from openbabel import *


def main(f, settings):

    """
    Find the axis atoms
    Find all the atoms to be rotated

    Rotate it and the substituents to the other side of the plane
    """
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()

    obconversion.ReadFile(obmol, f)

    obmol.ConnectTheDots()

    #Find the atoms composing furan ring
    Rings = obmol.GetSSSR()
    furan = []
    for ring in Rings:
        if len(settings.RingAtoms) == 5:
            if all(x in ring._path for x in settings.RingAtoms):
                furan = ring
                break
        else:
            if ring.Size() == 5 and not ring.IsAromatic():
                furan = ring
                break

    if furan == []:
        "No five membered rings to rotate. Quitting..."
        quit()
    #Find the plane of the 5-membered ring and the outlying atom
    norm, d, outAtom = FindFuranPlane(obmol, furan)

    #Find the atoms connected to the outlying atom and sort them
    #as either part of the ring(axis atoms) or as atoms to be rotated
    AxisAtoms = []
    RotAtoms = []

    for NbrAtom in OBAtomAtomIter(outAtom):
        #if NbrAtom.IsInRingSize(5):
        if furan.IsInRing(NbrAtom.GetIdx()):
            AxisAtoms.append(NbrAtom)
        else:
            RotAtoms.append(NbrAtom)
            FindSubstAtoms(NbrAtom, outAtom, RotAtoms)

    #Simple switch to help detect if the atoms are rotated the right way
    WasAbove90 = False
    angle = FindRotAngle(AxisAtoms[0], AxisAtoms[1], outAtom, norm)
    if angle > 0.5*pi:
        WasAbove90 = True
        rotangle = 2*(angle-0.5*pi)
    else:
        WasAbove90 = False
        rotangle = 2*(0.5*pi-angle)
    OldAtomCoords = outAtom.GetVector()
    print "Atom " + str(outAtom.GetAtomicNum()) + " will be rotated by " +\
        str(rotangle*57.3) + ' degrees'
    RotateAtom(outAtom, AxisAtoms[0], AxisAtoms[1], rotangle)
    angle2 = FindRotAngle(AxisAtoms[0], AxisAtoms[1], outAtom, norm)

    #if the atom is on the same side of the plane as it was,
    # it has been rotated in the wrong direction
    if ((angle2 > 0.5*pi) and WasAbove90) or ((angle2 < 0.5*pi) and not WasAbove90):
        #Flip the sign of the rotation angle, restore the coords
        #and rotate the atom in the opposite direction
        print "Atom was rotated the wrong way, switching the direction"
        rotangle = -rotangle
        outAtom.SetVector(OldAtomCoords)
        RotateAtom(outAtom, AxisAtoms[0], AxisAtoms[1], rotangle)

    RotatedAtoms = []  # Index to make sure that atoms are not rotated twice
    for atom in RotAtoms:
        if atom not in RotatedAtoms:
            RotateAtom(atom, AxisAtoms[0], AxisAtoms[1], rotangle)
            RotatedAtoms.append(atom)
        else:
            print "Atom already rotated, skipping"

    obconversion.SetOutFormat("sdf")
    obconversion.WriteFile(obmol, f[:-4] + 'rot.sdf')


#Recursively finds all the atoms connected to the input
def FindSubstAtoms(atom, outAtom, al):

    indexes = [a.GetIdx() for a in al]
    for NbrAtom in OBAtomAtomIter(atom):

        if (NbrAtom.GetIdx() not in indexes) and\
                (NbrAtom.GetIdx() != outAtom.GetIdx()):
            al.append(NbrAtom)
            FindSubstAtoms(NbrAtom, outAtom, al)


#Rotate atom around and axis by an angle
def RotateAtom(atom, AxisAtom1, AxisAtom2, angle):

    [u, v, w] = GetUnitVector(AxisAtom1, AxisAtom2)
    [x, y, z] = [atom.x(), atom.y(), atom.z()]
    [a, b, c] = [(AxisAtom1.x()+AxisAtom1.x())/2, (AxisAtom1.y()+AxisAtom1.y())/2,\
                 (AxisAtom1.z()+AxisAtom1.z())/2]

    X = (a*(v**2 + w**2) - u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(angle))+x*cos(angle)\
        +(-1*c*v+b*w-w*y+v*z)*sin(angle)
    Y = (b*(u**2 + w**2) - v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(angle))+y*cos(angle)\
        +(c*u-a*w+w*x-u*z)*sin(angle) #was _+_u*z)*sin(angle)
    Z = (c*(u**2 + v**2) - w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(angle))+z*cos(angle)\
        +(-1*b*u+a*v-v*x+u*y)*sin(angle)
    
    atom.SetVector(X, Y, Z)


def GetUnitVector(Atom1, Atom2):
    vector = []
    vector.append(Atom2.x() - Atom1.x())
    vector.append(Atom2.y() - Atom1.y())
    vector.append(Atom2.z() - Atom1.z())

    length = np.linalg.norm(vector)
    return [x/length for x in vector]


#Finds the angle by which atoms need to be rotated by taking the angle
#the atom is out of the plane (the 2 neighbor atoms being the axis)
#and doubling it
def FindRotAngle(AxisAtom1, AxisAtom2, OutAtom, Normal):

    start = []
    start.append((AxisAtom1.x() + AxisAtom2.x())/2)
    start.append((AxisAtom1.y() + AxisAtom2.y())/2)
    start.append((AxisAtom1.z() + AxisAtom2.z())/2)

    vector = []
    vector.append(OutAtom.x() - start[0])
    vector.append(OutAtom.y() - start[1])
    vector.append(OutAtom.z() - start[2])

    #Angle between plane normal and OOP atom
    vangle = angle(vector, Normal)

    #print "Measured angle: " + str(vangle*57.3)
    return vangle


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


"""Finds planes for every 3 atoms, calculates distances to the plane
for the other 2 atoms andchoose the plane with the smallest smallest distance
"""
def FindFuranPlane(mol, furan):

    atomIds = furan._path
    atoms = []

    for i in atomIds:
        atoms.append(mol.GetAtom(i))

    MinError = 100.0

    for atom in atoms:
        pats = [a for a in atoms if a != atom]
        norm, d, error = LstSqPlane(pats[0], pats[1], pats[2], pats[3])
        if error < MinError:
            MinError = error
            MaxNorm = norm
            MaxD = d
            OutAtom = atom

    return MaxNorm, MaxD, OutAtom


#Given 3 atoms, finds a plane defined by a normal vector and d
def FindPlane(atom1, atom2, atom3):

    vector1 = [atom2.x() - atom1.x(), atom2.y() - atom1.y(),
               atom2.z() - atom1.z()]
    vector2 = [atom3.x() - atom1.x(), atom3.y() - atom1.y(),
               atom3.z() - atom1.z()]
    cross_product = [vector1[1] * vector2[2] - vector1[2] * vector2[1],
                     -1 * vector1[0] * vector2[2] + vector1[2] * vector2[0],
                     vector1[0] * vector2[1] - vector1[1] * vector2[0]]

    d = cross_product[0] * atom1.x() - cross_product[1] * atom1.y() + \
        cross_product[2] * atom1.z()

    return cross_product, d


def LstSqPlane(atom1, atom2, atom3, atom4):

    # Inital guess of the plane
    [a0, b0, c0], d0 = FindPlane(atom1, atom2, atom3)

    f = lambda (a, b, c, d): PlaneError([atom1, atom2, atom3, atom4], a, b, c, d)
    res = sciopt.minimize(f, (a0, b0, c0, d0), method='nelder-mead')
    plane = list(res.x)

    return plane[:3], plane[3], f(plane)


def PlaneError(atoms, a, b, c, d):
    dists = []
    for atom in atoms:
        dists.append(abs(PointPlaneDist([a, b, c], d, atom)))
    return sum(dists)/len(dists)


#Calculates distance from an atom to a plane
def PointPlaneDist(norm, d, atom):

    point = []

    point.append(atom.x())
    point.append(atom.y())
    point.append(atom.z())

    a = norm[0]*point[0] + norm[1]*point[1] + norm[2]*point[2] + d
    b = sqrt(norm[0]**2 + norm[1]**2 + norm[2]**2)

    return a/b
