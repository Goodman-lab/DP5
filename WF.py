import numpy as np
import qml
from qml.fchl import get_atomic_kernels
import argparse
from scipy.stats import gaussian_kde as kde
import pickle
from scipy.stats import gmean
from pathlib import Path
from scipy import stats
import os
import pathos.multiprocessing as mp
import copy
from matplotlib import pyplot as plt

try:
    from openbabel.openbabel import OBConversion, OBMol, OBAtomAtomIter, OBMolAtomIter
except ImportError:
    from openbabel import *


c_distance = 4.532297920317418

class WFdata:

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

        self.UnscaledAtomProbs = [] #per atom unscaled wf probabilities for all conformers
        self.ScaledAtomProbs = [] #per atom wf scaled probabilities for all conformers

        self.CUnscaledprobs = []  #WF for isomers based on Carbon data
        self.CScaledprobs = []  # WF for isomers based on Carbon data
        self.Cplusprobs = []  # WF for isomers based on Carbon data

        self.BUnscaledAtomProbs = [] #per atom wf unscaled probabilities boltzmann weighted
        self.BScaledAtomProbs = [] #per atom wf scaled probabilities boltzmann weighted

        self.WFunscaledprobs = []  # Final WFS
        self.WFscaledprobs = []  # Final WFS
        self.WFplusprobs = []  # Final WFs plus probs

        self.output = str()  # final DP4 output

        if self.Atom_number < 86:

            self.folded_scaled_errors = pickle.load(open(ScriptPath /  "folded_scaled_errors.p", "rb"))

            self.atomic_reps = pickle.load(open(ScriptPath /  "atomic_reps_c.p", "rb"))

        else:

            self.folded_scaled_errors = pickle.load(open(ScriptPath / "frag_folded_errors.p", "rb"))

            self.atomic_reps = pickle.load(open(ScriptPath / "frag_reps.p", "rb"))

        self.mean_abs_error = np.mean(abs(self.folded_scaled_errors))

        self.output = ""


def kde_probs(Isomers,wfData,sigma):

    def kde_probfunction(conf_shifts, conf_reps):

        probs = []
        scaled_probs = []

        errors = [abs(shift - exp) for shift, exp in zip(conf_shifts, wfData.Cexp[iso])]

        scaled_shifts = ScaleNMR(conf_shifts, wfData.Cexp[iso])

        scaled_errors = [abs(shift - exp) for shift, exp in zip(scaled_shifts, wfData.Cexp[iso])]

        for e, s_e, r in zip(errors, scaled_errors, conf_reps):

            # calculate similarites between this atom and those in the atomic representation test set

            K_sim = get_atomic_kernels(np.array([r]), wfData.atomic_reps, [sigma],
                                       cut_distance=c_distance)[0][0]

            K_sim = np.hstack((K_sim, K_sim))

            # calculate kde using K_sim as the weighting function


            if np.sum(K_sim) == 0:

                kde_estimator = kde(wfData.folded_scaled_errors)

            else:

                kde_estimator = kde(wfData.folded_scaled_errors, weights=K_sim)

            e_diff = abs(e - wfData.mean_abs_error)

            p = kde_estimator.integrate_box_1d(wfData.mean_abs_error - e_diff, wfData.mean_abs_error + e_diff)

            probs.append(p)

            s_e_diff = abs(s_e - wfData.mean_abs_error)

            s_p = kde_estimator.integrate_box_1d(wfData.mean_abs_error - s_e_diff, wfData.mean_abs_error + s_e_diff)

            scaled_probs.append(s_p)

        return probs, scaled_probs

    #for each atom in the molecule calculate the atomic worry factor

    wfData.UnscaledAtomProbs = [[] for i in range(len(Isomers))]
    wfData.ScaledAtomProbs = [[] for i in range(len(Isomers))]

    for iso in range(len(Isomers)):

        res = [[] for i in wfData.AtomReps[iso]]
        wfData.UnscaledAtomProbs[iso] = [[] for i in wfData.AtomReps[iso]]
        wfData.ScaledAtomProbs[iso] = [[] for i in wfData.AtomReps[iso]]

        maxproc = 5

        pool = mp.Pool(maxproc)

        ind1 = 0

        for conf_shifts , conf_reps in zip(wfData.ConfCshifts[iso],wfData.AtomReps[iso] ) :
            res[ind1] = pool.apply_async(kde_probfunction,
                                         [conf_shifts,conf_reps])

            ind1 += 1

        for ind1 in range(len(res)):

            wfData.UnscaledAtomProbs[iso][ind1], wfData.ScaledAtomProbs[iso][ind1] = res[ind1].get()

    return wfData


def ProcessIsomers(WFdata, Isomers,Settings):


    OutputFolder = Path(Settings.OutputFolder)

    # extract calculated and experimental shifts and add to WFdata instance

    # Carbon

    # make sure any shifts with missing peaks are removed from all isomers

    removedC = []

    for iso in Isomers:

        WFdata.Cexp.append([])
        WFdata.Cshifts.append([])
        WFdata.Clabels.append([])

        WFdata.ConfCshifts.append([[] for i in range(len(iso.DFTConformers))])

        j = 0

        for shift, exp, label in zip(iso.Cshifts, iso.Cexp, iso.Clabels):

            if exp != '':

                WFdata.Cshifts[-1].append(shift)
                WFdata.Cexp[-1].append(exp)
                WFdata.Clabels[-1].append(label)

                for i in range( len(WFdata.ConfCshifts[-1])):

                    WFdata.ConfCshifts[-1][i].append(iso.ConformerCShifts[i][j])

                    i+=1

            elif label not in removedC:

                removedC.append(label)

            j+=1

    for l in removedC:

        for j, Clabel in enumerate(WFdata.Clabels):

            if l in Clabel:
                i = Clabel.index(l)

                WFdata.Cshifts[j].pop(i)

                WFdata.Cexp[j].pop(i)

                WFdata.Clabels[j].pop(i)


    #write qml compound objects and atomic representations



    #check the number of atoms in the structures

    #if there are less than 86 (max number of atoms in a molecule in the training set) atoms

    if WFdata.Atom_number < 86:

        for iso in Isomers:

            #open new xyz file

            InputFile = Path(iso.InputFile)

            #find conformer with the lowest energy

            WFdata.AtomReps.append([])

            for i,geom in enumerate(iso.DFTConformers):

                xyz_file = open(str(OutputFolder / "wf" /InputFile.stem) + "_" +str(i).zfill(3) + ".xyz", "w")

                xyz_file.write(str(len(iso.Atoms)) + "\n" + "\n")

                for atom, coords in zip(iso.Atoms, geom):

                    xyz_file.write(atom + " " + str(coords[0]) + " " + str(coords[1]) + " " + str(coords[2]) + "\n")

                xyz_file.close()

                WFdata.Compounds.append(qml.Compound(xyz = str(Settings.OutputFolder/"wf"/ InputFile.stem) +"_"+ str(i).zfill(3) + ".xyz"))

                WFdata.Compounds[-1].generate_fchl_representation(max_size=86, cut_distance=c_distance)

                WFdata.AtomReps[-1].append([])

                for C_l in iso.Clabels:

                    ind = int(C_l.split("C")[1])

                    WFdata.AtomReps[-1][-1].append(WFdata.Compounds[-1].representation[ind])

    #otherwise we need to fragment the molecule to radius of 3

    else:

        for iso in Isomers:

            #open new xyz file

            InputFile = Path(iso.InputFile)

            #find conformer with the lowest energy

            WFdata.AtomReps.append([])

            for i,geom in enumerate(iso.DFTConformers):

                xyz_file = open(str(OutputFolder / "wf" /InputFile.stem) + "_" +str(i).zfill(3) + ".xyz", "w")

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

                os.mkdir(str(OutputFolder / "wf" /InputFile.stem) + "_" +str(i).zfill(3) + "_fragments")

                mol_fragments(m,str(OutputFolder / "wf" /InputFile.stem) + "_" +str(i).zfill(3) + "_fragments")

                conf_rep = []

                for xyz_frag in sorted(os.listdir(  str(OutputFolder / "wf" /InputFile.stem) + "_" +str(i).zfill(3) + "_fragments")):

                    c = qml.Compound(xyz=str(OutputFolder / "wf" /InputFile.stem) + "_" +str(i).zfill(3) + "_fragments/" + xyz_frag)

                    c.generate_fchl_representation(max_size=54, cut_distance=c_distance)

                    conf_rep.append(c.representation[0])

                WFdata.AtomReps[-1].append([])

                for C_l in iso.Clabels:

                    ind = int(C_l.split("C")[1])

                    WFdata.AtomReps[-1][-1].append(conf_rep[ind])

    return WFdata


def InternalScaling(WFdata):
    # perform internal scaling process

    # calculate prediction errors

    if len(WFdata.Cexp[0]) > 0:

        for Cshifts, Cexp in zip(WFdata.Cshifts, WFdata.Cexp):
            WFdata.Cscaled.append(ScaleNMR(Cshifts, Cexp))

    '''

    if len(WFdata.Hexp[0]) > 0:

        for Hshifts, Hexp in zip(WFdata.Hshifts, WFdata.Hexp):
            WFdata.Hscaled.append(ScaleNMR(Hshifts, Hexp))

        for Hscaled, Hexp in zip(WFdata.Hscaled, WFdata.Hexp):
            WFdata.Hscalederrors.append([Hscaled[i] - Hexp[i] for i in range(0, len(Hscaled))])
    '''
    return WFdata


def ScaleNMR(calcShifts, expShifts):

    slope, intercept, r_value, p_value, std_err = stats.linregress(expShifts,
                                                                   calcShifts)
    scaled = [(x - intercept) / slope for x in calcShifts]

    return scaled


def BoltzmannWeight_WF(Isomers,WFdata):

    print("Conf", np.shape(np.array(WFdata.ScaledAtomProbs)))

    for iso,scaled_probs,unscaled_probs in zip( Isomers, WFdata.ScaledAtomProbs,WFdata.UnscaledAtomProbs):

        B_scaled_probs = [0] * len(scaled_probs[0])

        B_unscaled_probs = [0] * len(scaled_probs[0])

        for population, conf_scaled_p,conf_unscaled_p in zip(iso.Populations, scaled_probs,unscaled_probs ):

            for i in range(len(B_scaled_probs)):

                B_scaled_probs[i] += conf_scaled_p[i] * population

                B_unscaled_probs[i] += conf_unscaled_p[i] * population

        WFdata.BScaledAtomProbs.append(B_scaled_probs)

        WFdata.BUnscaledAtomProbs.append(B_unscaled_probs)

    return WFdata


def Calculate_WF(WFdata):

    for scaled_probs,unscaled_probs in zip(WFdata.BScaledAtomProbs,WFdata.BUnscaledAtomProbs):

        WFunscaled = 1 - gmean([1 - p_si for p_si in unscaled_probs])

        WFscaled = (1 - gmean([1 - p_si for p_si in scaled_probs]))

        WFplus = (1  - gmean([1 - p_si for p_si in scaled_probs] + [1 - p_si for p_si in unscaled_probs]))

        WFdata.CUnscaledprobs.append(WFunscaled)

        WFdata.CScaledprobs.append(WFscaled)

        WFdata.Cplusprobs.append(WFplus)

    return WFdata


def Rescale_WF(WFdata,Settings):

    incorrect_kde = pickle.load(open(Path(Settings.ScriptDir) / "i_w_kde_mean_s_0.025.p" ,"rb"))

    correct_kde = pickle.load(open(Path(Settings.ScriptDir) / "c_w_kde_mean_s_0.025.p" ,"rb"))

    i = 0

    for scaled,unscaled in zip(WFdata.BScaledAtomProbs,WFdata.BUnscaledAtomProbs):

        WFdata.BScaledAtomProbs[i] = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in scaled  ]
        WFdata.BUnscaledAtomProbs[i] = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in unscaled  ]

        i += 1

    WFdata.WFunscaledprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in WFdata.CUnscaledprobs  ]  # Final WFS
    WFdata.WFscaledprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in WFdata.CScaledprobs  ]   # Final WFS
    WFdata.WFplusprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in WFdata.Cplusprobs  ]   # Final WFs plus probs

    return WFdata


def Pickle_res(WFdata,Settings):

    data_dic = {"Cshifts": WFdata.Cshifts,
                "Cexp": WFdata.Cexp,
                "Clabels": WFdata.Clabels,
                "Hshifts": WFdata.Hshifts,
                "Hexp": WFdata.Hexp,
                "Hlabels": WFdata.Hlabels,
                "Cscaled": WFdata.Cscaled,
                "Hscaled": WFdata.Hscaled,
                "ConfCshifts": WFdata.ConfCshifts,
                "Compounds": WFdata.Compounds,
                "AtomReps": WFdata.AtomReps,
                "UnscaledAtomProbs": WFdata.UnscaledAtomProbs,
                "ScaledAtomProbs": WFdata.ScaledAtomProbs,
                "BUnscaledAtomProbs": WFdata.BUnscaledAtomProbs,
                "BScaledAtomProbs": WFdata.BScaledAtomProbs,
                "CUnscaledprobs": WFdata.CUnscaledprobs,
                "CScaledprobs": WFdata.CScaledprobs,
                "Cplusprobs": WFdata.Cplusprobs,
                "WFunscaledprobs": WFdata.WFunscaledprobs,
                "WFscaledprobs": WFdata.WFscaledprobs,
                "WFplusprobs": WFdata.WFplusprobs}

    pickle.dump(data_dic , open(Path(Settings.OutputFolder) / "wf" / "data_dic.p","wb"))

    return WFdata


def UnPickle_res(WFdata,Settings):

    data_dic =  pickle.load(open(Path(Settings.OutputFolder) / "wf" / "data_dic.p","rb"))

    WFdata.Cshifts = data_dic["Cshifts"]
    WFdata.Cexp = data_dic["Cexp"]
    WFdata.Clabels = data_dic["Clabels"]
    WFdata.Hshifts = data_dic["Hshifts"]
    WFdata.Hexp = data_dic["Hexp"]
    WFdata.Hlabels = data_dic["Hlabels"]
    WFdata.Cscaled = data_dic["Cscaled"]
    WFdata.Hscaled = data_dic["Hscaled"]
    WFdata.ConfCshifts = data_dic["ConfCshifts"]
    WFdata.Compounds = data_dic["Compounds"]
    WFdata.AtomReps = data_dic["AtomReps"]
    WFdata.UnscaledAtomProbs = data_dic["UnscaledAtomProbs"]
    WFdata.ScaledAtomProbs = data_dic["ScaledAtomProbs"]
    WFdata.BUnscaledAtomProbs = data_dic["BUnscaledAtomProbs"]
    WFdata.BScaledAtomProbs = data_dic["BScaledAtomProbs"]
    WFdata.CUnscaledprobs = data_dic["CUnscaledprobs"]
    WFdata.CScaledprobs = data_dic["CScaledprobs"]
    WFdata.Cplusprobs = data_dic["Cplusprobs"]
    WFdata.WFunscaledprobs = data_dic["WFunscaledprobs"]
    WFdata.WFscaledprobs = data_dic["WFscaledprobs"]
    WFdata.WFplusprobs = data_dic["WFplusprobs"]

    return WFdata


def PrintAssignment(WFData):
    isomer = 0

    for Clabels, Cshifts, Cexp, Cscaled, atom_p in zip(WFData.Clabels, WFData.Cshifts, WFData.Cexp, WFData.Cscaled,WFData.BScaledAtomProbs):
        WFData.output += ("\n\nAssigned C shifts for isomer " + str(isomer + 1) + ": ")

        PrintNMR(Clabels, Cshifts, Cscaled, Cexp,atom_p, WFData)

        isomer += 1


def PrintNMR(labels, values, scaled, exp,atom_p, WFData):
    s = np.argsort(values)

    svalues = np.array(values)[s]

    slabels = np.array(labels)[s]

    sscaled = np.array(scaled)[s]

    sexp = np.array(exp)[s]

    atom_p = np.array(atom_p)[s]

    WFData.output += ("\nlabel, calc, corrected, exp, error,prob")

    for i in range(len(labels)):
        WFData.output += ("\n" + format(slabels[i], "6s") + ' ' + format(svalues[i], "6.2f") + ' '
                           + format(sscaled[i], "6.2f") + ' ' + format(sexp[i], "6.2f") + ' ' +
                           format(sexp[i] - sscaled[i], "6.2f")+ ' ' +
                           format(atom_p[i], "6.2f"))


def MakeOutput(WFData, Isomers, Settings):
    # add some info about the calculation

    WFData.output += Settings.InputFiles[0] + "\n"

    WFData.output += "\n" + "Solvent = " + Settings.Solvent

    WFData.output += "\n" + "Force Field = " + Settings.ForceField + "\n"

    if 'o' in Settings.Workflow:
        WFData.output += "\n" + "DFT optimisation Functional = " + Settings.oFunctional
        WFData.output += "\n" + "DFT optimisation Basis = " + Settings.oBasisSet

    if 'e' in Settings.Workflow:
        WFData.output += "\n" + "DFT energy Functional = " + Settings.eFunctional
        WFData.output += "\n" + "DFT energy Basis = " + Settings.eBasisSet

    if 'n' in Settings.Workflow:
        WFData.output += "\n" + "DFT NMR Functional = " + Settings.nFunctional
        WFData.output += "\n" + "DFT NMR Basis = " + Settings.nBasisSet

    if Settings.StatsParamFile != "none":
        WFData.output += "\n\nStats model = " + Settings.StatsParamFile

    WFData.output += "\n\nNumber of isomers = " + str(len(Isomers))

    c = 1

    for i in Isomers:
        WFData.output += "\nNumber of conformers for isomer " + str(c) + " = " + str(len(i.Conformers))

        c += 1

    PrintAssignment(WFData)

    WFData.output += ("\n\nResults of DP4 using Carbon: ")

    for i, p in enumerate(WFData.WFscaledprobs):
        WFData.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    WFData.output += ("\n\nResults of DP4: ")


    print("number of c carbons = " + str(len(Isomers[0].Clabels)))
    print("number of e carbons = " + str(len(WFData.Cexp[0])))

    print(WFData.output)

    if Settings.OutputFolder == '':

        out = open(str(os.getcwd()) + "/" + str(Settings.InputFiles[0] + "NMR.wf"), "w+")

    else:

        out = open(os.path.join(Settings.OutputFolder, str(Settings.InputFiles[0] + "NMR.wf")), "w+")

    out.write(WFData.output)

    out.close()

    return WFData


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