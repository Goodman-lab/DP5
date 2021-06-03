# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 15:20:18 2014

@author: ke291

Code for diastereomer, tautomer and protomer generation via InChI strings.
This file gets called by PyDP4.py if diastereomer and/or tautomer and/or
protomer generation is used.
"""
from PyDP4 import settings
import sys
import os

try:
    from openbabel.openbabel import OBConversion, OBMol, OBAtomAtomIter, OBMolAtomIter
except ImportError:
    from openbabel import *

import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

def main(f):

    inchi, aux = GetInchi(f)

    ds_inchis = GenDiastereomers(inchi)
    ds_inchis = [FixTautProtons(f, i, aux) for i in ds_inchis]

    for ds in range(0, len(ds_inchis)):

        print("Isomer " + str(ds) + " inchi = " + ds_inchis[ds])

        Inchi2Struct(ds_inchis[ds], f[:-4] + str(ds+1), aux)
        RestoreNumsSDF(f[:-4] + str(ds+1) + '.sdf', f, aux)


def GetInchiRenumMap(AuxInfo):

    for l in AuxInfo.split('/'):
        if 'N:' in l:
            RenumLayer = l
            break
    amap = [int(x) for x in RenumLayer[2:].split(',')]
    return amap


def FixTautProtons(f, inchi, AuxInfo):

    #Get tautomeric protons and atoms they are connected to from Inchi
    TautProts = GetTautProtons(inchi)
    amap = GetInchiRenumMap(AuxInfo)

    #get the correspondence of the Inchi numbers to the source numbers
    hmap = []
    for taut in TautProts:
        for heavyatom in range(1, len(taut)):
            hmap.append([int(taut[heavyatom]), amap[int(taut[heavyatom])-1]])

    #Read molecule from file
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, f)

    Fixprotpos = []
    for heavyatom in hmap:
        atom = obmol.GetAtom(heavyatom[1])
        for nbratom in OBAtomAtomIter(atom):
            if nbratom.GetAtomicNum() == 1:
                Fixprotpos.append(heavyatom[0])
    draftFH = []
    for i in range(0, len(Fixprotpos)):
        if Fixprotpos[i] not in [a[0] for a in draftFH]:
            draftFH.append([Fixprotpos[i], Fixprotpos.count(Fixprotpos[i])])

    fixedlayer = '/f/h'
    for h in draftFH:
        if h[1] == 1:
            fixedlayer = fixedlayer + str(h[0])+'H,'
        else:
            fixedlayer = fixedlayer + str(h[0])+'H' + str(h[1]) + ','

    resinchi = inchi + fixedlayer[:-1]

    return resinchi


#Get H connections from sdf file
def GetHcons(f):
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, f)
    Hcons = []
    for atom in OBMolAtomIter(obmol):
        idx = atom.GetIdx()
        anum = atom.GetAtomicNum()
        if anum == 1:
            for NbrAtom in OBAtomAtomIter(atom):
                Hcons.append([idx, NbrAtom.GetIdx()])
    return Hcons


def RestoreNumsSDF(f, fold, AuxInfo):

    #Read molecule from file
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, f)
    #Get the atoms Hs are connected to
    oldHcons = GetHcons(fold)
    #translate the H connected atoms to the new numbering system
    amap = GetInchiRenumMap(AuxInfo)
    for i in range(0, len(oldHcons)):
        oldHcons[i][1] = amap.index(oldHcons[i][1])+1

    newHcons = []
    temp = []
    i = 0
    for atom in OBMolAtomIter(obmol):
        idx = atom.GetIdx()
        anum = atom.GetAtomicNum()
        #If atom is hydrogen, check what it is connected to
        if anum == 1:
            for NbrAtom in OBAtomAtomIter(atom):
                newHcons.append([idx, NbrAtom.GetIdx()])
        #Pick the temporary atom
        temp.append(atom)

    for i in range(0, len(newHcons)):
        conatom = newHcons[i][1]
        for b in range(0, len(oldHcons)):
            if conatom == oldHcons[b][1]:
                amap.append(oldHcons[b][0])
                #remove the number, so that it doesn't get added twice
                oldHcons[b][1] = 0

    newmol = OBMol()
    added = []

    for i in range(1, len(amap)+1):
        newn = amap.index(i)
        newmol.AddAtom(temp[newn])
        added.append(newn)

    #Final runthrough to check that all atoms have been added,
    #tautomeric protons can be missed. If tautomeric proton tracking
    #is implemented this can be removed
    for i in range(0, len(temp)):
        if not i in added:
            newmol.AddAtom(temp[i])

    #Restore the bonds
    newmol.ConnectTheDots()
    newmol.PerceiveBondOrders()
    #Write renumbered molecule to file
    obconversion.SetOutFormat("sdf")
    obconversion.WriteFile(newmol, f)


def GetInchi(f):

    print("Getting inchi from file ",f)

    if os.path.sep not in f:
        f = os.path.join(os.getcwd(), f)

    m = Chem.MolFromMolFile(f, removeHs = False)

    m = Chem.AddHs(m)

    idata = Chem.MolToInchiAndAuxInfo(m)

    return idata[0], idata[1]


def Inchi2Struct(inchi, f, aux):


    cwd = os.getcwd()
    fullf = os.path.join(cwd, f)
    infile = open(f + '.inchi', 'w')
    infile.write(inchi)
    infile.close()

    inchi = open(f + '.inchi', "r").read()

    m = AllChem.inchi.MolFromInchi(inchi, sanitize=True, removeHs=False)

    m = AllChem.AddHs(m, addCoords=True)

    AllChem.EmbedMolecule(m)

    save3d = Chem.SDWriter(fullf + '.sdf')

    save3d.write(m)


def GetTautProtons(inchi):
    #get the tautomer layer and pickup the data
    layers = inchi.split('/')

    for l in layers:
        if 'h' in l:
            ProtLayer = l
    ProtList = list(ProtLayer)
    starts = []
    ends = []
    for i in range(0, len(ProtList)):
        if ProtList[i] == '(':
            starts.append(i)
        if ProtList[i] == ')':
            ends.append(i)
    TautProts = []
    for i in range(0, len(starts)):
        TautProts.append((ProtLayer[starts[i]+1:ends[i]]).split(','))

    return TautProts


def GenSelectDiastereomers(structf, atoms):

    f = structf

    if (f[-4:] != '.sdf'):
        f += '.sdf'

    inchi, aux = GetInchi(f)
    amap = GetInchiRenumMap(aux)

    translated_atoms = []
    for atom in atoms:
        translated_atoms.append(amap.index(atom)+1)

    ds_inchis = GenSelectDSInchis(inchi, translated_atoms)
    ds_inchis = [FixTautProtons(f, i, aux) for i in ds_inchis]
    filenames = []
    for ds in range(0, len(ds_inchis)):
        Inchi2Struct(ds_inchis[ds], f[:-4] + str(ds + 1), aux)
        RestoreNumsSDF(f[:-4] + str(ds + 1) + '.sdf', f, aux)
        filenames.append(f[:-4] + str(ds + 1))

    return filenames


def GenSelectDSInchis(inchi, atoms):
    #Inchis of all diastereomers, including the parent structure
    resinchis = []

    #get the number of potential diastereomers
    layers = inchi.decode().split('/')
    for l in layers:
        if 't' in l:
            slayer = l
            sc = l[1:].decode().split(',')

    ignore = []
    for i in range(0, len(sc)):
        if not int(sc[i][:-1]) in atoms:
            ignore.append(sc[i])
    sc = [x for x in sc if x not in ignore]

    if len(sc) == 0:
        "No stereocentres remaining, no diastereomers will be generated."
        return 0

    numds = 2**(len(sc))
    print("Number of diastereomers to be generated: " + str(numds))
    temps = []
    #Generate inversion patterns - essentially just binary strings
    for i in range(0, numds):
        template = bin(i)[2:].zfill(len(sc))
        temps.append(template)

    #For each 1 in template, invert the corresponding stereocentre
    #and add the resulting diastereomer to the list
    invert = {'+': '-', '-': '+'}

    reslayers = []
    for ds in range(0, numds):
        newds = list(sc)
        for stereocentre in range(0, len(sc)):
            if temps[ds][stereocentre] == '1':
                tlist = list(newds[stereocentre])
                tlist[-1] = invert[tlist[-1]]
                newds[stereocentre] = "".join(tlist)
        newlayer = str(slayer)
        for stereocentre in range(0, len(sc)):
            newlayer = newlayer.replace(sc[stereocentre], newds[stereocentre])
        reslayers.append(newlayer)
    print(reslayers)
    resinchis = []
    for layer in reslayers:
        resinchis.append(inchi.replace(slayer, layer))
    return resinchis


def GenDiastereomers(structf, nS, atoms=[]):
    if len(atoms) > 0:
        return GenSelectDiastereomers(structf, atoms)

    f = structf

    if (f[-4:] != '.sdf'):
        f += '.sdf'

    if nS < 2:
        cwd = os.getcwd()

        fullf = os.path.join(cwd, f)

        shutil.copy(fullf, fullf[:-4] + "0.sdf")

        return [f[:-4] + "0.sdf"]

    inchi, aux = GetInchi(f)


    i,a = GetInchi(f)

    ds_inchis = GenDSInchis(inchi)
    ds_inchis = [FixTautProtons(f, i, aux) for i in ds_inchis]
    filenames = []

    for ds in range(0, len(ds_inchis)):

        print("Isomer " + str(ds) + " inchi = " + ds_inchis[ds])

        Inchi2Struct(ds_inchis[ds], f[:-4] + str(ds + 1), aux)
        RestoreNumsSDF(f[:-4] + str(ds + 1) + '.sdf', f, aux)
        filenames.append(f[:-4] + str(ds + 1))

    return filenames


def GenDSInchis(inchi):

    ilist = list(inchi)
    #Inchis of all diastereomers, including the parent structure
    resinchis = []

    #get the number of potential diastereomers
    numds = 0
    layers = inchi.split('/')
    for l in layers:
        if 't' in l:
            numds = 2**(len(l.translate({ord(i): None for i in 't,1234567890'}))-1)

    if numds == 0:

        raise ValueError("No chiral carbon detected in the input molecule!")

    else:
        print("Number of diastereomers to be generated: " + str(numds))

    #find configuration sites (+ and -)
    bs = ilist.index('t')
    es = ilist[bs:].index('/')
    spos = []
    for s in range(bs, bs+es):
        if ilist[s] == '+' or ilist[s] == '-':
            spos.append(s)

    temps = []
    #Generate inversion patterns - essentially just binary strings
    for i in range(0, numds):
        template = bin(i)[2:].zfill(len(spos)-1)
        temps.append(template)

    #For each 1 in template, invert the corresponding stereocentre
    #and add the resulting diastereomer to the list
    invert = {'+': '-', '-': '+'}

    for ds in range(0, numds):
        t = list(ilist)
        for stereocentre in range(1, len(spos)):
            if temps[ds][stereocentre-1] == '1':
                t[spos[stereocentre]] = invert[t[spos[stereocentre]]]
        resinchis.append(''.join(t))

    return resinchis
