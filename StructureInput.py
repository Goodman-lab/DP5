from rdkit import Chem
from rdkit.Chem import AllChem
import os
from InchiGen import *

def GenerateMolFromSDF(InputFile):

    cwd = os.getcwd()

    fullf = os.path.join(cwd, InputFile + '.sdf')

    m = Chem.MolFromMolFile(fullf, removeHs=False,sanitize = True)

    return m

def GenerateSDFFromTxt(InputFile,inp_type):

    Mols = []

    cwd = os.getcwd()

    InputFile

    if os.path.exists(os.path.join(cwd, InputFile)):

        fullf = os.path.join(cwd, InputFile)

    elif inp_type == 'InChI':

        fullf = os.path.join(cwd, InputFile + ".inchi")

    elif inp_type == 'Smiles':

        fullf = os.path.join(cwd, InputFile + ".smiles")

    elif inp_type == 'Smarts':

        fullf = os.path.join(cwd, InputFile + ".smarts")

    f = open(fullf, "r")

    if inp_type == 'InChI':

        for line in f.readlines():

            if line.strip() == '':
                continue

            else:

                try:

                    m = Chem.MolFromInchi(line.strip(),sanitize = True)

                    Mols.append(m)

                except:

                    print(line + " could not be read")

    elif inp_type == 'Smiles':

        for line in f.readlines():

            if line.strip() == '':
                continue

            else:

                try:

                    m = Chem.MolFromSmiles(line.strip(),sanitize = True)

                    Mols.append(m)

                except:

                    print(line + " could not be read")

    elif inp_type == 'Smarts':

        for line in f.readlines():

            if line.strip() == '':
                continue

            else:

                try:

                    m = Chem.MolFromSmarts(line.strip(),sanitize = True)

                    Mols.append(m)

                except:

                    print(line + " could not be read")

    else:

        print('unrecognised')

    GeneratedFiles = GenerateSDFFromMols(Mols,inp_type )

    return GeneratedFiles

def GenerateSDFFromMols(mols,inp_type):

    files = []

    cwd = os.getcwd()

    for ind, m in enumerate(mols):

        m = AllChem.AddHs(m, addCoords=True)

        AllChem.EmbedMolecule(m)

        AllChem.MMFFOptimizeMolecule(m)

        Chem.rdmolops.AssignStereochemistryFrom3D(m)

        f = inp_type +  '_Mol_' + str(ind) +"_.sdf"

        files.append(f[:-4])

        fullf = os.path.join(cwd, f)

        save3d = Chem.SDWriter( fullf)

        save3d.write(m)

    return files

def CleanUp(InputFiles):

    # check input file types

    CleanedInputFiles = []

    cwd = os.getcwd()

    for f in InputFiles:

        if f.endswith('.sdf'):

            f = f[:-4]

        fullf = os.path.join(cwd, f + 'cleaned.sdf')

        m = GenerateMolFromSDF(f)

        m = AllChem.AddHs(m, addCoords=True)

        AllChem.EmbedMolecule(m)

        AllChem.MMFFOptimizeMolecule(m)

        Chem.rdmolops.AssignStereochemistryFrom3D(m)

        save3d = Chem.SDWriter(fullf)

        save3d.write(m)

        CleanedInputFiles.append( f + 'cleaned')

    return CleanedInputFiles

def NumberofStereoCentres(InputFile):

    cwd = os.getcwd()

    fullf = os.path.join(cwd, InputFile + '.sdf')

    m = Chem.MolFromMolFile(fullf, removeHs=False)

    m = AllChem.AddHs(m, addCoords=True)

    AllChem.EmbedMolecule(m)

    AllChem.MMFFOptimizeMolecule(m)

    Chem.rdmolops.AssignStereochemistryFrom3D(m)

    nStereo = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(m)

    return nStereo