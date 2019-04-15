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
sys.path.append('/home/' + settings.user + '/Tools/openbabel-install/lib/python2.7/site-packages/')
import os
from openbabel import *
import subprocess
import itertools

MolConPath = '/home/' + settings.user + '/chemaxon/marvinsuite/bin/molconvert'


def main(f):

    inchi, aux = GetInchi(f)

    print(inchi)

    ds_inchis = GenDiastereomers(inchi)
    ds_inchis = [FixTautProtons(f, i, aux) for i in ds_inchis]

    for ds in range(0, len(ds_inchis)):
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
    
    cwd = os.getcwd()
    outp = subprocess.check_output(MolConPath + ' inchi ' + cwd + '/' + f, shell=True)
    idata = outp.split('\n')

    aux = idata[1][:]
    return idata[0], aux


def Inchi2Struct(inchi, f, aux):

    cwd = os.getcwd()
    fullf = cwd + '/' + f
    infile = open(f + '.inchi', 'w')
    infile.write(inchi)
    infile.close()
    
    outp = subprocess.check_output(MolConPath + ' sdf ' + fullf +
        '.inchi -3:c1S{fine}[prehydrogenize] -o ' + fullf + '.sdf', shell=True)


def GenProtomers(structf, atoms):

    f = structf + ".sdf"
    inchi, aux = GetInchi(f)
    print(inchi)
    amap = GetInchiRenumMap(aux)

    prot_atoms = []
    for atom in atoms:
        prot_atoms.append(amap.index(atom))

    finchi = FixTautProtons(f, inchi, aux)
    prot_inchis = GenProtInchis(finchi, prot_atoms)

    print(prot_inchis)
    filenames = []
    for prot in range(0, len(prot_inchis)):
        Inchi2Struct(prot_inchis[prot], f[:-4] + 'p' + str(prot+1), aux)
        RestoreNumsSDF(f[:-4] + 'p' + str(prot+1) + '.sdf', f, aux)
        filenames.append(f[:-4] + 'p' + str(prot+1))
    return len(filenames), filenames


def GenProtInchis(inchi, atoms):

    #Read and interpret proton counts on heavy atoms from inchi,
    # including tautomeric layers
    protons, formula, tautomerics, fprotons = ReadProtonCounts(inchi)

    #Construct list of heavy atoms with tautomeric protons and which system
    #they belong to
    tlist = [[], []]
    i = 0
    for tlayer in tautomerics:
        temp = tlayer[1:-1].split(',')
        tlist[0].extend([int(x) for x in temp[1:]])
        tlist[1].extend([i for x in temp[1:]])
    print(tlist)
    #Increase H count and regenerate the formula
    for i in range(0, len(formula)):
        if formula[i][0] == 'H':
            formula[i][1] += 1
        formula[i][1] = str(formula[i][1])
    formula = [''.join(x) for x in formula]
    formula = ''.join(formula)

    #For each basic atom in atoms, generate a copy of original protons
    #add atom and save it
    print("Protonating atoms with these InChI numbers: " +\
        str([x+1 for x in atoms]))
    protlayers = []
    fprotlayers = []
    for atom in atoms:
        if atom+1 not in tlist[0]:
            extraprotons = list(protons)
            extrafprotons = list(fprotons)
            extraprotons[atom] += 1
            if tautomerics == []:
                protlayers.append(WriteProtonCounts(extraprotons))
            else:
                protlayers.append(WriteProtonCounts(extraprotons) +
                                  ',' + ','.join(tautomerics))
            fprotlayers.append(WriteProtonCounts(extrafprotons))
        else:
            extraprotons = list(protons)
            extrafprotons = list(fprotons)
            extratautomerics = list(tautomerics)
            extrafprotons[atom] += 1
            #which tautomeric system atom belongs to?
            tindex = tlist[0].index(atom+1)
            tindex = tlist[1][tindex]
            temp = tautomerics[tindex].split(',')
            #Get the proton count and increase by 1
            protcount = int(temp[0][2:]) + 1
            #Write the proton count back
            temp[0] = temp[0][:2] + str(protcount)
            extratautomerics[tindex] = ','.join(temp)
            protlayers.append(WriteProtonCounts(extraprotons)+',' +
                              ','.join(extratautomerics))
            fprotlayers.append(WriteProtonCounts(extrafprotons))
    protinchis = []
    protinchis.append(inchi)
    for l in range(0, len(protlayers)):
        layers = inchi.split('/')
        MainLayerPassed = False
        ChargeAdded = False
        i = 1
        while (i < len(layers)):
            if 'h' in layers[i]:
                if not MainLayerPassed:
                    layers[i] = protlayers[l]
                    MainLayerPassed = True
                    if 'q' not in inchi:
                        layers.insert(i+1, 'q+1')
                        ChargeAdded = True
                else:
                    layers[i] = fprotlayers[l]
            if ('q' in layers[i]) and ChargeAdded is False:
                charge = int(layers[i][1:])
                layers[i] = 'q'+"%+d" % charge
            if 'C' in layers[i] and 'H' in layers[i]:
                layers[i] = formula
            i += 1
        #insert charge layer here
        protinchis.append('/'.join(layers))
    return protinchis


def WriteProtonCounts(protons):
    collectedprotons = [[], [], [], []]

    i = 1
    lastcount = protons[0]
    start = 0
    while i < len(protons):
        if protons[i] != lastcount:
            if start == i-1:
                collectedprotons[lastcount].append(str(start+1))
            else:
                collectedprotons[lastcount].append(str(start+1)+'-'+str(i))
            lastcount = protons[i]
            start = i
        i += 1

    if start == i-1:
        collectedprotons[lastcount].append(str(start+1))
    else:
        collectedprotons[lastcount].append(str(start+1)+'-'+str(i))

    hlayer = 'h'
    if len(collectedprotons[1]) > 0:
        hlayer += ','.join(collectedprotons[1])
        hlayer += 'H,'
    if len(collectedprotons[2]) > 0:
        hlayer += ','.join(collectedprotons[2])
        hlayer += 'H2,'
    if len(collectedprotons[3]) > 0:
        hlayer += ','.join(collectedprotons[3])
        hlayer += 'H3,'
    hlayer = hlayer[:-1]
    return hlayer


def ReadProtonCounts(inchi):
    import re

    #Get inchi layers
    layers = inchi.split('/')
    ProtLayer = ''
    FixedLayer = ''
    for l in layers[1:]:
        if 'C' in l and 'H' in l:
            atoms = re.findall(r"[a-zA-Z]+", l)
            indexes = [int(x) for x in re.findall(r"\d+", l)]
            formula = [list(x) for x in zip(atoms, indexes)]
        if 'h' in l and ProtLayer != '':
            FixedLayer = l[1:]
        if 'h' in l and ProtLayer == '':
            ProtLayer = l[1:]

    #initialize proton list
    nheavy = sum([x[1] for x in formula if x[0] != 'H'])

    #Find, save and remove tautomeric portions from main proton layer
    tautomerics = re.findall(r"\(.*?\)", ProtLayer)
    ProtLayer = re.sub(r"\(.*?\)", "", ProtLayer)
    if ProtLayer[-1] == ',':
        ProtLayer = ProtLayer[:-1]

    #Read the main and the fixed proton layer
    protons = ReadPSections(ProtLayer, nheavy)
    fprotons = ReadPSections(FixedLayer, nheavy)

    return protons, formula, tautomerics, fprotons


def ReadPSections(ProtLayer, nheavy):
    import re
    protons = [0 for x in range(0, nheavy)]
    #seperate the 1proton, 2proton and 3proton atoms, then seperate the records
    psections = [x for x in re.findall(r".*?(?=H)", ProtLayer) if x != '']
    secvals = [0 for x in range(0, len(psections))]
    psections[0] = psections[0].split(',')

    #interpret each record and fill in the proton table
    #start by finding the proton count value for each section
    for i in range(1, len(psections)):
        if psections[i][0] == ',':
            secvals[i-1] = 1
            psections[i] = psections[i][1:].split(',')
        else:
            secvals[i-1] = int(psections[i][0])
            psections[i] = psections[i][2:].split(',')
    if ProtLayer[-1] != 'H':
        secvals[-1] = int(ProtLayer[-1])
    else:
        secvals[-1] = 1

    #now expand each entry in the sections and fill the corresponding value
    #in proton table
    for i in range(0, len(psections)):
        for s in psections[i]:
            if '-' in s:
                [start, finish] = [int(x) for x in s.split('-')]
                protons[start-1:finish] = [secvals[i] for x in
                    range(0, len(protons[start-1:finish]))]
            else:
                protons[int(s)-1] = secvals[i]
    return protons


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


def GenTautomers(structf):

    f = structf + ".sdf"
    inchi, aux = GetInchi(f)

    abc = 'abcdefghijklmnopqrstuvwxyz'

    t_inchis = GenTautInchis(inchi, aux, structf)
    filenames = []
    for ds in range(0, len(t_inchis)):
        Inchi2Struct(t_inchis[ds], f[:-4] + abc[ds], aux)
        RestoreNumsSDF(f[:-4] + abc[ds] + '.sdf', f, aux)
        filenames.append(f[:-4] + abc[ds])
    return len(filenames), filenames


#For now only works on one tautomeric system - for nucleosides don't need more
def GenTautInchis(inchi, aux, structf):

    resinchis = []  # Inchis of all tautomers, including the parent structure

    #Get tautomeric protons and atoms they are connected to
    TautProts = GetTautProtons(inchi)

    #The total number of tautomeric protons in the system
    if str(TautProts[0][0][1:]) != '':
        totprotons = int(str(TautProts[0][0][1:]))
    else:
        totprotons = 1

    #Check for and remove non-hetero atoms
    temp = [int(x) for x in TautProts[0][1:] if IsHetero(int(x), inchi)]

    #Get the numbering map, to see which atomes are the tautomeric ones in the
    #original structure file. From there determine their type and valency
    #based on connectivity
    amap = GetInchiRenumMap(aux)
    OldTautProts = [amap[prot-1] for prot in temp]
    valencies = GetTautValency(structf, OldTautProts)

    #the multivalent atoms will always have at least one proton
    superfixed = []
    for i in range(0, len(valencies)):
        if valencies[i] > 1:
            superfixed.append(TautProts[0][i+1])
            #TautProts.append(['H', TautProts[0][i+1]])

    #Generate all the possible proton positions
    #with repetitions for multivalency
    fixedprotons = list(itertools.combinations(TautProts[0][1:],
                        r=totprotons-len(superfixed)))

    for i in range(0, len(fixedprotons)):
        fixedprotons[i] = superfixed + list(fixedprotons[i])
        fixedprotons[i].sort(key=int)

    #Count the occurences of protons positions, save the counts and
    #remove redundant positions
    counts = []
    for i in range(0, len(fixedprotons)):
        counts.append([])
        j = 0
        while j < len(fixedprotons[i]):
            counts[i].append(fixedprotons[i].count(fixedprotons[i][j]))
            if counts[i][-1] > 1:
                for _ in range(0, counts[i][-1]-1):
                    fixedprotons[i].remove(fixedprotons[i][j])
            #j+=counts[i][-1]
            j += 1

    fixprefix = '/f/h'
    tauts = []
    for i in range(0, len(fixedprotons)):
        tauts.append(fixprefix)
        for j in range(0, len(fixedprotons[i])):
            if j > 0:
                tauts[i] += ','
            tauts[i] += fixedprotons[i][j] + 'H'
            if counts[i][j] > 1:
                tauts[i] += str(counts[i][j])

    for taut in tauts:
        resinchis.append(inchi + taut)

    return resinchis


def GetTautValency(structf, tautatoms):
     #Read molecule from file
    obconversion = OBConversion()
    obconversion.SetInFormat("sdf")
    obmol = OBMol()
    obconversion.ReadFile(obmol, structf + ".sdf")

    protvalency = []
    for idx in tautatoms:
        atom = obmol.GetAtom(idx)
        corenum = atom.GetAtomicNum()
        #If atom is hydrogen, check what it is connected to
        nheavynbr = 0
        for NbrAtom in OBAtomAtomIter(atom):
            nbrnum = NbrAtom.GetAtomicNum()
            if nbrnum > 1:
                nheavynbr += 1
        if corenum == 8:
            protvalency.append(1)
        if corenum == 7:
            protvalency.append(3-nheavynbr)

    return protvalency


#utility function that determines if an atom is a heteroatom
def IsHetero(n, inchi):
    layers = inchi.split('/')
    nC = int((layers[1].split('H'))[0][1:])
    if (n > nC):
        return True
    else:
        return False


def GenSelectDiastereomers(structf, atoms):

    f = structf + ".sdf"
    inchi, aux = GetInchi(f)
    amap = GetInchiRenumMap(aux)

    translated_atoms = []
    for atom in atoms:
        translated_atoms.append(amap.index(atom)+1)

    ds_inchis = GenSelectDSInchis(inchi, translated_atoms)
    ds_inchis = [FixTautProtons(f, i, aux) for i in ds_inchis]
    filenames = []
    for ds in range(0, len(ds_inchis)):
        Inchi2Struct(ds_inchis[ds], f[:-4] + str(ds+1), aux)
        RestoreNumsSDF(f[:-4] + str(ds+1) + '.sdf', f, aux)
        filenames.append(f[:-4] + str(ds+1))
    return len(filenames), filenames


def GenSelectDSInchis(inchi, atoms):
    #Inchis of all diastereomers, including the parent structure
    resinchis = []

    #get the number of potential diastereomers
    layers = inchi.split('/')
    for l in layers:
        if 't' in l:
            slayer = l
            sc = l[1:].split(',')

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


def GenDiastereomers(structf):

    f = structf + ".sdf"
    inchi, aux = GetInchi(f)

    print(inchi)

    ds_inchis = GenDSInchis(inchi)
    ds_inchis = [FixTautProtons(f, i, aux) for i in ds_inchis]
    filenames = []
    for ds in range(0, len(ds_inchis)):
        Inchi2Struct(ds_inchis[ds], f[:-4] + str(ds+1), aux)
        RestoreNumsSDF(f[:-4] + str(ds+1) + '.sdf', f, aux)
        filenames.append(f[:-4] + str(ds+1))
    return len(filenames), filenames


def GenDSInchis(inchi):

    ilist = list(inchi)
    #Inchis of all diastereomers, including the parent structure
    resinchis = []

    #get the number of potential diastereomers
    layers = inchi.split('/')
    for l in layers:
        if 't' in l:
            numds = 2**(len(l.translate(None, 't,1234567890'))-1)

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
