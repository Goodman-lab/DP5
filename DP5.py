import numpy as np
import qml
from qml.fchl import get_atomic_kernels
from scipy.stats import gaussian_kde as kde
import pickle
from scipy.stats import gmean
from pathlib import Path
from scipy import stats
import os
import pathos.multiprocessing as mp
import copy
import gzip

try:
    from openbabel.openbabel import OBConversion, OBMol, OBAtomAtomIter, OBMolAtomIter
except ImportError:
    from openbabel import *


c_distance = 4.532297920317418

class DP5data:

    def __init__(self,ScriptPath,Atoms):

        self.Atom_number = Atoms

        self.Cshifts = []  # Carbon shifts used in DP4 calculation
        self.Cexp = []  # Carbon experimental shifts used in DP4 calculation
        self.Clabels = []  # Carbon atom labels
        self.Hshifts = []  # Proton shifts used in DP4 calculation
        self.Hexp = []  # Proton experimental shifts used in DP4 calculation
        self.Hlabels = []  # Proton atom labels
        self.Cscaled = []  # Internally scaled carbon shifts
        self.Hscaled = []  # Internally scaled proton shifts

        self.ConfCshifts = []

        self.Compounds = [] #qml compound objects for each isomer
        self.AtomReps = [] # FCHL representations ordered by Clabels

        self.ScaledAtomProbs = [] #per atom dp5 scaled probabilities for all conformers

        self.CScaledprobs = []  # DP5 for isomers based on Carbon data

        self.BScaledAtomProbs = [] #per atom dp5 scaled probabilities boltzmann weighted

        self.DP5scaledprobs = []  # Final DP5S

        self.output = str()  # final DP4 output

        self.folded_scaled_errors = pickle.load(open(ScriptPath / "folded_scaled_errors.p", "rb"))

        #self.folded_unscaled_errors = pickle.load(open(ScriptPath / "folded_unscaled_errors.p", "rb"))

        if self.Atom_number < 86:

            with gzip.open(ScriptPath / "atomic_reps_c.gz", "rb") as f:

                self.atomic_reps = pickle.load(f)

        else:

            with gzip.open(ScriptPath / "frag_reps.gz", "rb") as f:

                self.atomic_reps = pickle.load(f)

        self.mean_abs_error = np.mean(abs(self.folded_scaled_errors))

        self.output = ""

def kde_probs(Isomers,dp5Data,sigma):

    def kde_probfunction(conf_shifts, conf_reps):

        scaled_probs = []

        errors = [abs(shift - exp) for shift, exp in zip(conf_shifts, dp5Data.Cexp[iso])]

        scaled_shifts = ScaleNMR(conf_shifts, dp5Data.Cexp[iso])

        scaled_errors = [abs(shift - exp) for shift, exp in zip(scaled_shifts, dp5Data.Cexp[iso])]

        for e, s_e, r in zip(errors, scaled_errors, conf_reps):

            # calculate similarites between this atom and those in the atomic representation test set

            K_sim = get_atomic_kernels(np.array([r]), dp5Data.atomic_reps, [sigma],
                                       cut_distance=c_distance)[0][0]

            K_sim = np.hstack((K_sim, K_sim))

            # calculate kde using K_sim as the weighting function

            if np.sum(K_sim) == 0:

                scaled_kde_estimator = kde(dp5Data.folded_scaled_errors)

            else:

                scaled_kde_estimator = kde(dp5Data.folded_scaled_errors, weights=K_sim)

            s_e_diff = abs(s_e - dp5Data.mean_abs_error)

            s_p = scaled_kde_estimator.integrate_box_1d(dp5Data.mean_abs_error - s_e_diff, dp5Data.mean_abs_error + s_e_diff)

            scaled_probs.append(s_p)

        return  scaled_probs

    #for each atom in the molecule calculate the atomic worry factor

    dp5Data.ScaledAtomProbs = [[] for i in range(len(Isomers))]

    for iso in range(len(Isomers)):

        res = [[] for i in dp5Data.AtomReps[iso]]

        dp5Data.ScaledAtomProbs[iso] = [[] for i in dp5Data.AtomReps[iso]]

        maxproc = 4

        pool = mp.Pool(maxproc)

        ind1 = 0

        for conf_shifts , conf_reps in zip(dp5Data.ConfCshifts[iso],dp5Data.AtomReps[iso] ) :
            res[ind1] = pool.apply_async(kde_probfunction,
                                         [conf_shifts,conf_reps])

            ind1 += 1

        for ind1 in range(len(res)):

            dp5Data.ScaledAtomProbs[iso][ind1] = res[ind1].get()

    return dp5Data


def ProcessIsomers(dp5Data, Isomers,Settings):


    OutputFolder = Path(Settings.OutputFolder)

    # extract calculated and experimental shifts and add to dp5Data instance

    # Carbon

    # make sure any shifts with missing peaks are removed from all isomers

    removedC = []

    for iso in Isomers:

        dp5Data.Cexp.append([])
        dp5Data.Cshifts.append([])
        dp5Data.Clabels.append([])

        dp5Data.ConfCshifts.append([[] for i in range(len(iso.DFTConformers))])

        j = 0

        for shift, exp, label in zip(iso.Cshifts, iso.Cexp, iso.Clabels):

            if exp != '':

                dp5Data.Cshifts[-1].append(shift)
                dp5Data.Cexp[-1].append(exp)
                dp5Data.Clabels[-1].append(label)

                for i in range( len(dp5Data.ConfCshifts[-1])):

                    dp5Data.ConfCshifts[-1][i].append(iso.ConformerCShifts[i][j])

                    i+=1

            elif label not in removedC:

                removedC.append(label)

            j+=1

    for l in removedC:

        for j, Clabel in enumerate(dp5Data.Clabels):

            if l in Clabel:
                i = Clabel.index(l)

                dp5Data.Cshifts[j].pop(i)

                dp5Data.Cexp[j].pop(i)

                dp5Data.Clabels[j].pop(i)


    #write qml compound objects and atomic representations

    #check the number of atoms in the structures

    #if there are less than 86 (max number of atoms in a molecule in the training set) atoms

    if dp5Data.Atom_number < 86:

        for iso in Isomers:

            #open new xyz file

            InputFile = Path(iso.InputFile)

            #find conformer with the lowest energy

            dp5Data.AtomReps.append([])

            for i,geom in enumerate(iso.DFTConformers):

                xyz_file = open(str(OutputFolder / "dp5" /InputFile.stem) + "_" +str(i).zfill(3) + ".xyz", "w")

                xyz_file.write(str(len(iso.Atoms)) + "\n" + "\n")

                for atom, coords in zip(iso.Atoms, geom):

                    xyz_file.write(atom + " " + str(coords[0]) + " " + str(coords[1]) + " " + str(coords[2]) + "\n")

                xyz_file.close()

                dp5Data.Compounds.append(qml.Compound(xyz = str(Settings.OutputFolder/"dp5"/ InputFile.stem) +"_"+ str(i).zfill(3) + ".xyz"))

                dp5Data.Compounds[-1].generate_fchl_representation(max_size=86, cut_distance=c_distance)

                dp5Data.AtomReps[-1].append([])

                for C_l in iso.Clabels:

                    ind = int(C_l.split("C")[1])

                    dp5Data.AtomReps[-1][-1].append(dp5Data.Compounds[-1].representation[ind])

    #otherwise we need to fragment the molecule to radius of 3

    else:

        for iso in Isomers:

            #open new xyz file

            InputFile = Path(iso.InputFile)

            #find conformer with the lowest energy

            dp5Data.AtomReps.append([])

            for i,geom in enumerate(iso.DFTConformers):

                xyz_file = open(str(OutputFolder / "dp5" /InputFile.stem) + "_" +str(i).zfill(3) + ".xyz", "w")

                xyz_file.write(str(len(iso.Atoms)) + "\n" + "\n")

                for atom, coords in zip(iso.Atoms, geom):

                    xyz_file.write(atom + " " + str(coords[0]) + " " + str(coords[1]) + " " + str(coords[2]) + "\n")

                xyz_file.close()

                #now need to fragment the molecule and generate these representations

                #build ob mol

                obconversion = OBConversion()
                obconversion.SetInFormat("sdf")
                m = OBMol()
                obconversion.ReadFile(m, iso.InputFile)

                os.mkdir(str(OutputFolder / "dp5" /InputFile.stem) + "_" +str(i).zfill(3) + "_fragments")

                mol_fragments(m,str(OutputFolder / "dp5" /InputFile.stem) + "_" +str(i).zfill(3) + "_fragments")

                conf_rep = []

                for xyz_frag in sorted(os.listdir(  str(OutputFolder / "dp5" /InputFile.stem) + "_" +str(i).zfill(3) + "_fragments")):

                    c = qml.Compound(xyz=str(OutputFolder / "dp5" /InputFile.stem) + "_" +str(i).zfill(3) + "_fragments/" + xyz_frag)

                    c.generate_fchl_representation(max_size=54, cut_distance=c_distance)

                    conf_rep.append(c.representation[0])

                dp5Data.AtomReps[-1].append([])

                for C_l in iso.Clabels:

                    ind = int(C_l.split("C")[1])

                    dp5Data.AtomReps[-1][-1].append(conf_rep[ind])

    return dp5Data


def InternalScaling(dp5Data):
    # perform internal scaling process

    # calculate prediction errors

    if len(dp5Data.Cexp[0]) > 0:

        for Cshifts, Cexp in zip(dp5Data.Cshifts, dp5Data.Cexp):
            dp5Data.Cscaled.append(ScaleNMR(Cshifts, Cexp))

    '''

    if len(dp5Data.Hexp[0]) > 0:

        for Hshifts, Hexp in zip(dp5Data.Hshifts, dp5Data.Hexp):
            dp5Data.Hscaled.append(ScaleNMR(Hshifts, Hexp))

        for Hscaled, Hexp in zip(dp5Data.Hscaled, dp5Data.Hexp):
            dp5Data.Hscalederrors.append([Hscaled[i] - Hexp[i] for i in range(0, len(Hscaled))])
    '''
    return dp5Data


def ScaleNMR(calcShifts, expShifts):

    slope, intercept, r_value, p_value, std_err = stats.linregress(expShifts,
                                                                   calcShifts)
    scaled = [(x - intercept) / slope for x in calcShifts]

    return scaled


def BoltzmannWeight_DP5(Isomers,dp5Data):

    for iso,scaled_probs in zip( Isomers, dp5Data.ScaledAtomProbs):

        B_scaled_probs = [0] * len(scaled_probs[0])

        for population, conf_scaled_p in zip(iso.Populations, scaled_probs ):

            for i in range(len(B_scaled_probs)):

                B_scaled_probs[i] += conf_scaled_p[i] * population

        dp5Data.BScaledAtomProbs.append(B_scaled_probs)

    return dp5Data


def Calculate_DP5(dp5Data):

    for scaled_probs in dp5Data.BScaledAtomProbs:

        DP5scaled = 1 - gmean([1 - p_si for p_si in scaled_probs])

        dp5Data.CScaledprobs.append(DP5scaled)

    return dp5Data


def Rescale_DP5(dp5Data,Settings):

    incorrect_kde = pickle.load(open(Path(Settings.ScriptDir) / "i_w_kde_mean_s_0.025.p" ,"rb"))

    correct_kde = pickle.load(open(Path(Settings.ScriptDir) / "c_w_kde_mean_s_0.025.p" ,"rb"))

    i = 0

    for scaled in dp5Data.BScaledAtomProbs:

        dp5Data.BScaledAtomProbs[i] = [ 1 - float(incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x))) for x in scaled  ]

        i += 1

    dp5Data.DP5scaledprobs = [  1 - float(incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x))) for x in dp5Data.CScaledprobs  ]   # Final DP5S

    return dp5Data


def Pickle_res(dp5Data,Settings):

    data_dic = {"Cshifts": dp5Data.Cshifts,
                "Cexp": dp5Data.Cexp,
                "Clabels": dp5Data.Clabels,
                "Hshifts": dp5Data.Hshifts,
                "Hexp": dp5Data.Hexp,
                "Hlabels": dp5Data.Hlabels,
                "Cscaled": dp5Data.Cscaled,
                "Hscaled": dp5Data.Hscaled,
                "ConfCshifts": dp5Data.ConfCshifts,
                "Compounds": dp5Data.Compounds,
                "AtomReps": dp5Data.AtomReps,

                "ScaledAtomProbs": dp5Data.ScaledAtomProbs,

                "BScaledAtomProbs": dp5Data.BScaledAtomProbs,

                "CScaledprobs": dp5Data.CScaledprobs,

                "DP5scaledprobs": dp5Data.DP5scaledprobs}

    pickle.dump(data_dic , open(Path(Settings.OutputFolder) / "dp5" / "data_dic.p","wb"))

    return dp5Data


def UnPickle_res(dp5Data,Settings):

    data_dic =  pickle.load(open(Path(Settings.OutputFolder) / "dp5" / "data_dic.p","rb"))

    dp5Data.Cshifts = data_dic["Cshifts"]
    dp5Data.Cexp = data_dic["Cexp"]
    dp5Data.Clabels = data_dic["Clabels"]
    dp5Data.Hshifts = data_dic["Hshifts"]
    dp5Data.Hexp = data_dic["Hexp"]
    dp5Data.Hlabels = data_dic["Hlabels"]
    dp5Data.Cscaled = data_dic["Cscaled"]
    dp5Data.Hscaled = data_dic["Hscaled"]
    dp5Data.ConfCshifts = data_dic["ConfCshifts"]
    dp5Data.Compounds = data_dic["Compounds"]
    dp5Data.AtomReps = data_dic["AtomReps"]
    dp5Data.ScaledAtomProbs = data_dic["ScaledAtomProbs"]
    dp5Data.BScaledAtomProbs = data_dic["BScaledAtomProbs"]
    dp5Data.CScaledprobs = data_dic["CScaledprobs"]
    dp5Data.DP5scaledprobs = data_dic["DP5scaledprobs"]

    return dp5Data


def PrintAssignment(dp5Data):
    isomer = 0

    for Clabels, Cshifts, Cexp, Cscaled, atom_p in zip(dp5Data.Clabels, dp5Data.Cshifts, dp5Data.Cexp, dp5Data.Cscaled,dp5Data.BScaledAtomProbs):
        dp5Data.output += ("\n\nAssigned C shifts for isomer " + str(isomer + 1) + ": ")

        PrintNMR(Clabels, Cshifts, Cscaled, Cexp,atom_p, dp5Data)

        isomer += 1


def PrintNMR(labels, values, scaled, exp,atom_p, dp5Data):

    s = np.argsort(values)

    svalues = np.array(values)[s]

    slabels = np.array(labels)[s]

    sscaled = np.array(scaled)[s]

    sexp = np.array(exp)[s]

    atom_p = np.array(atom_p)[s]

    dp5Data.output += ("\nlabel, calc, corrected, exp, error,prob")

    for i in range(len(labels)):

        dp5Data.output += ("\n" + format(slabels[i], "6s") + ' ' + format(svalues[i], "6.2f") + ' '
                           + format(sscaled[i], "6.2f") + ' ' + format(sexp[i], "6.2f") + ' ' +
                           format(sexp[i] - sscaled[i], "6.2f")+ ' ' +
                           format(atom_p[i] , "6.2f"))


def MakeOutput(dp5Data, Isomers, Settings):
    # add some info about the calculation

    dp5Data.output += Settings.InputFiles[0] + "\n"

    dp5Data.output += "\n" + "Solvent = " + Settings.Solvent

    dp5Data.output += "\n" + "Force Field = " + Settings.ForceField + "\n"

    if 'o' in Settings.Workflow:
        dp5Data.output += "\n" + "DFT optimisation Functional = " + Settings.oFunctional
        dp5Data.output += "\n" + "DFT optimisation Basis = " + Settings.oBasisSet

    if 'e' in Settings.Workflow:
        dp5Data.output += "\n" + "DFT energy Functional = " + Settings.eFunctional
        dp5Data.output += "\n" + "DFT energy Basis = " + Settings.eBasisSet

    if 'n' in Settings.Workflow:
        dp5Data.output += "\n" + "DFT NMR Functional = " + Settings.nFunctional
        dp5Data.output += "\n" + "DFT NMR Basis = " + Settings.nBasisSet

    if Settings.StatsParamFile != "none":
        dp5Data.output += "\n\nStats model = " + Settings.StatsParamFile

    dp5Data.output += "\n\nNumber of isomers = " + str(len(Isomers))

    c = 1

    for i in Isomers:
        dp5Data.output += "\nNumber of conformers for isomer " + str(c) + " = " + str(len(i.Conformers))

        c += 1

    PrintAssignment(dp5Data)

    dp5Data.output += ("\n\nResults of DP5 using Carbon: ")

    for i, p in enumerate(dp5Data.DP5scaledprobs):

        dp5Data.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    print(dp5Data.output)

    if Settings.OutputFolder == '':

        out = open(str(os.getcwd()) + "/" + str(Settings.InputFiles[0] + "NMR.dp5"), "w+")

    else:

        out = open(os.path.join(Settings.OutputFolder, str(Settings.InputFiles[0] + "NMR.dp5")), "w+")

    out.write(dp5Data.output)

    out.close()

    return dp5Data


def mol_fragments(mole,outfile):

    obconv = OBConversion()

    obconv.SetOutFormat("xyz")

    c = 1

    for atom, exp,dft,diff in zip(OBMolAtomIter(mole)):

        # if this is a carbon atom start a breadth first search for other carbon atoms with depth specified

        # create a new mol instance

        new_mol = OBMol()

        # add this atom

        # new_mol.AddAtom(atom)

        fragment_ind = []

        l = atom.GetIndex()

        fragment_ind.append(l)

        # for iteration depth radius

        old_queue = [atom]

        for iteration in range(0, 3):

            new_queue = []

            for a in old_queue:

                for atom2 in OBAtomAtomIter(a):

                    i = atom2.GetIndex()

                    # if the atom has not been seen before add it to the fragment ind list and to the new molecule

                    if i not in fragment_ind:
                        new_queue.append(atom2)

                        fragment_ind.append(i)

                        # new_mol.AddAtom(atom2)

            old_queue = copy.copy(new_queue)

        fragment_ind = [fragment_ind[0]] + sorted(fragment_ind[1:])

        for i in fragment_ind:

            for a in OBMolAtomIter(mole):

                if a.GetIndex() == i:
                    new_mol.AddAtom(a)

        f = open(outfile + "frag" + str(l).zfill(3) + ".xyz", "w+")

        f.write(str(new_mol.NumAtoms()) + "\n\n")

        i = 0

        for atom in OBMolAtomIter(new_mol):

            f.write(atom.GetType()[0] + " " + str(atom.GetX()) + " " + str(atom.GetY()) + " " + str(
                atom.GetZ()) + "\n")


            i+=1

        f.close()

        c += 1