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
        self.ConfHshifts = []

        self.Compounds = [] #qml compound objects for each isomer
        self.CAtomReps = [] # FCHL representations ordered by Clabels
        self.HAtomReps = [] # FCHL representations ordered by Hlabels

        self.CUnscaledAtomProbs = [] #per atom unscaled wf probabilities for all conformers
        self.CScaledAtomProbs = [] #per atom wf scaled probabilities for all conformers

        self.HUnscaledAtomProbs = [] #per atom unscaled wf probabilities for all conformers
        self.HScaledAtomProbs = [] #per atom wf scaled probabilities for all conformers

        self.CUnscaledprobs = []  #DP5 for isomers based on Carbon data
        self.CScaledprobs = []  # DP5 for isomers based on Carbon data
        self.Cplusprobs = []  # DP5 for isomers based on Carbon data

        self.HUnscaledprobs = []  #DP5 for isomers based on Carbon data
        self.HScaledprobs = []  # DP5 for isomers based on Carbon data
        self.Hplusprobs = []  # DP5 for isomers based on Carbon data

        self.CBUnscaledAtomProbs = [] #per atom wf unscaled probabilities boltzmann weighted
        self.CBScaledAtomProbs = [] #per atom wf scaled probabilities boltzmann weighted

        self.HBUnscaledAtomProbs = [] #per atom wf unscaled probabilities boltzmann weighted
        self.HBScaledAtomProbs = [] #per atom wf scaled probabilities boltzmann weighted

        self.CDP5unscaledprobs = []  # Final DP5S
        self.CDP5scaledprobs = []  # Final DP5S
        self.CDP5plusprobs = []  # Final DP5s plus probs

        self.HDP5unscaledprobs = []  # Final DP5S
        self.HDP5scaledprobs = []  # Final DP5S
        self.HDP5plusprobs = []  # Final DP5s plus probs

        self.DP5unscaledprobs = []  # Final combined
        self.DP5scaledprobs = []  # Final combined
        self.DP5plusprobs = []  # Final combined

        self.output = str()  # final DP4 output

        if self.Atom_number < 86:

            self.Cfolded_scaled_errors = pickle.load(open(ScriptPath /  "Cfolded_scaled_errors.p", "rb"))
            self.Hfolded_scaled_errors = pickle.load(open(ScriptPath /  "Hfolded_scaled_errors.p", "rb"))
            with gzip.open(ScriptPath / "Hatomic_reps_c.gz", "rb") as f:

                self.Hatomic_reps = pickle.load(f)

            with gzip.open(ScriptPath / "Catomic_reps_c.gz", "rb") as f:

                self.Catomic_reps = pickle.load(f)

        else:

            self.Cfolded_scaled_errors = pickle.load(open(ScriptPath / "Cfrag_folded_errors.p", "rb"))
            self.Hfolded_scaled_errors = pickle.load(open(ScriptPath / "Hfrag_folded_errors.p", "rb"))
            with gzip.open(ScriptPath / "Cfrag_reps.gz", "rb") as f:

                self.Catomic_reps = pickle.load(f)

            with gzip.open(ScriptPath / "Hfrag_reps.gz", "rb") as f:

                self.Hatomic_reps = pickle.load(f)

        self.Hmean_abs_error = np.mean(abs(self.Hfolded_scaled_errors))
        self.Cmean_abs_error = np.mean(abs(self.Cfolded_scaled_errors))
        self.output = ""


def kde_probs(Isomers,dp5Data,sigma):

    def kde_probfunction(Cconf_shifts, Cconf_reps,Hconf_shifts, Hconf_reps):

        Cprobs = []
        Cscaled_probs = []

        Cerrors = [abs(shift - exp) for shift, exp in zip(Cconf_shifts, dp5Data.Cexp[iso])]

        Cscaled_shifts = ScaleNMR(Cconf_shifts, dp5Data.Cexp[iso])

        Cscaled_errors = [abs(shift - exp) for shift, exp in zip(Cscaled_shifts, dp5Data.Cexp[iso])]

        for e, s_e, r in zip(Cerrors, Cscaled_errors, Cconf_reps):

            # calculate similarites between this atom and those in the atomic representation test set

            K_sim = get_atomic_kernels(np.array([r]), dp5Data.Catomic_reps, [sigma],
                                       cut_distance=c_distance)[0][0]

            K_sim = np.hstack((K_sim, K_sim))

            # calculate kde using K_sim as the weighting function


            if np.sum(K_sim) == 0:

                kde_estimator = kde(dp5Data.Cfolded_scaled_errors)

            else:

                kde_estimator = kde(dp5Data.Cfolded_scaled_errors, weights=K_sim)

            e_diff = abs(e - dp5Data.Cmean_abs_error)

            p = kde_estimator.integrate_box_1d(dp5Data.Cmean_abs_error - e_diff, dp5Data.Cmean_abs_error + e_diff)

            Cprobs.append(p)

            s_e_diff = abs(s_e - dp5Data.Cmean_abs_error)

            s_p = kde_estimator.integrate_box_1d(dp5Data.Cmean_abs_error - s_e_diff, dp5Data.Cmean_abs_error + s_e_diff)

            Cscaled_probs.append(s_p)

        Hprobs = []

        Hscaled_probs = []

        Herrors = [abs(shift - exp) for shift, exp in zip(Hconf_shifts, dp5Data.Hexp[iso])]

        Hscaled_shifts = ScaleNMR(Hconf_shifts, dp5Data.Hexp[iso])

        Hscaled_errors = [abs(shift - exp) for shift, exp in zip(Hscaled_shifts, dp5Data.Hexp[iso])]

        for e, s_e, r in zip(Herrors, Hscaled_errors, Hconf_reps):

            # calculate similarites between this atom and those in the atomic representation test set

            K_sim = get_atomic_kernels(np.array([r]), dp5Data.Hatomic_reps, [sigma],
                                       cut_distance=c_distance)[0][0]

            K_sim = np.hstack((K_sim, K_sim))

            # calculate kde using K_sim as the weighting function


            if np.sum(K_sim) == 0:

                kde_estimator = kde(dp5Data.Hfolded_scaled_errors)

            else:

                kde_estimator = kde(dp5Data.Hfolded_scaled_errors, weights=K_sim)

            e_diff = abs(e - dp5Data.Hmean_abs_error)

            p = kde_estimator.integrate_box_1d(dp5Data.Hmean_abs_error - e_diff, dp5Data.Hmean_abs_error + e_diff)

            Hprobs.append(p)

            s_e_diff = abs(s_e - dp5Data.Hmean_abs_error)

            s_p = kde_estimator.integrate_box_1d(dp5Data.Hmean_abs_error - s_e_diff, dp5Data.Hmean_abs_error + s_e_diff)

            Hscaled_probs.append(s_p)

        return Cprobs, Cscaled_probs,Hprobs, Hscaled_probs

    #for each atom in the molecule calculate the atomic worry factor

    dp5Data.CUnscaledAtomProbs = [[] for i in range(len(Isomers))]
    dp5Data.CScaledAtomProbs = [[] for i in range(len(Isomers))]

    dp5Data.HUnscaledAtomProbs = [[] for i in range(len(Isomers))]
    dp5Data.HScaledAtomProbs = [[] for i in range(len(Isomers))]


    for iso in range(len(Isomers)):

        res = [[] for i in dp5Data.AtomReps[iso]]

        dp5Data.CUnscaledAtomProbs[iso] = [[] for i in dp5Data.CAtomReps[iso]]
        dp5Data.CScaledAtomProbs[iso] = [[] for i in dp5Data.CAtomReps[iso]]

        dp5Data.HUnscaledAtomProbs[iso] = [[] for i in dp5Data.HAtomReps[iso]]
        dp5Data.HScaledAtomProbs[iso] = [[] for i in dp5Data.HAtomReps[iso]]

        maxproc = 5

        pool = mp.Pool(maxproc)

        ind1 = 0

        for Cconf_shifts , Cconf_reps,Hconf_shifts , Hconf_reps in zip(dp5Data.ConfCshifts[iso],dp5Data.CAtomReps[iso],dp5Data.ConfHshifts[iso],dp5Data.HAtomReps[iso] ) :
            res[ind1] = pool.apply_async(kde_probfunction,
                                         [Cconf_shifts,Cconf_reps,Hconf_shifts,Hconf_reps])

            ind1 += 1

        for ind1 in range(len(res)):

            dp5Data.CUnscaledAtomProbs[iso][ind1], dp5Data.CScaledAtomProbs[iso][ind1],dp5Data.HUnscaledAtomProbs[iso][ind1], dp5Data.HScaledAtomProbs[iso][ind1] = res[ind1].get()

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

    removedH = []

    for iso in Isomers:

        dp5Data.Hexp.append([])
        dp5Data.Hshifts.append([])
        dp5Data.Hlabels.append([])

        dp5Data.ConfHshifts.append([[] for i in range(len(iso.DFTConformers))])

        j = 0

        for shift, exp, label in zip(iso.Hshifts, iso.Hexp, iso.Hlabels):

            if exp != '':

                dp5Data.Hshifts[-1].append(shift)
                dp5Data.Hexp[-1].append(exp)
                dp5Data.Hlabels[-1].append(label)

                for i in range(len(dp5Data.ConfHshifts[-1])):
                    dp5Data.ConfHshifts[-1][i].append(iso.ConformerHShifts[i][j])

                    i += 1

            elif label not in removedH:

                removedH.append(label)

            j += 1

    for l in removedH:

        for j, Hlabel in enumerate(dp5Data.Hlabels):

            if l in Hlabel:
                i = Hlabel.index(l)

                dp5Data.Hshifts[j].pop(i)

                dp5Data.Hexp[j].pop(i)

                dp5Data.Hlabels[j].pop(i)


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

                xyz_file = open(str(OutputFolder / "wf" /InputFile.stem) + "_" +str(i).zfill(3) + ".xyz", "w")

                xyz_file.write(str(len(iso.Atoms)) + "\n" + "\n")

                for atom, coords in zip(iso.Atoms, geom):

                    xyz_file.write(atom + " " + str(coords[0]) + " " + str(coords[1]) + " " + str(coords[2]) + "\n")

                xyz_file.close()

                dp5Data.Compounds.append(qml.Compound(xyz = str(Settings.OutputFolder/"wf"/ InputFile.stem) +"_"+ str(i).zfill(3) + ".xyz"))

                dp5Data.Compounds[-1].generate_fchl_representation(max_size=86, cut_distance=c_distance)

                dp5Data.CAtomReps[-1].append([])

                for C_l in iso.Clabels:

                    ind = int(C_l.split("C")[1])

                    dp5Data.CAtomReps[-1][-1].append(dp5Data.CCompounds[-1].representation[ind])

                dp5Data.HAtomReps[-1].append([])

                for H_l in iso.Hlabels:

                    ind = int(H_l.split("H")[1])

                    dp5Data.HAtomReps[-1][-1].append(dp5Data.Compounds[-1].representation[ind])

    #otherwise we need to fragment the molecule to radius of 3

    else:

        for iso in Isomers:

            #open new xyz file

            InputFile = Path(iso.InputFile)

            #find conformer with the lowest energy

            dp5Data.AtomReps.append([])

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

                dp5Data.CAtomReps[-1].append([])
                dp5Data.HAtomReps[-1].append([])

                for C_l in iso.Clabels:

                    ind = int(C_l.split("C")[1])

                    dp5Data.CAtomReps[-1][-1].append(conf_rep[ind])

                for H_l in iso.Hlabels:

                    ind = int(H_l.split("C")[1])

                    dp5Data.HAtomReps[-1][-1].append(conf_rep[ind])

    return dp5Data


def InternalScaling(dp5Data):
    # perform internal scaling process

    # calculate prediction errors

    if len(dp5Data.Cexp[0]) > 0:

        for Cshifts, Cexp in zip(dp5Data.Cshifts, dp5Data.Cexp):
            dp5Data.Cscaled.append(ScaleNMR(Cshifts, Cexp))


    if len(dp5Data.Hexp[0]) > 0:

        for Hshifts, Hexp in zip(dp5Data.Hshifts, dp5Data.Hexp):
            dp5Data.Hscaled.append(ScaleNMR(Hshifts, Hexp))

    return dp5Data


def ScaleNMR(calcShifts, expShifts):

    slope, intercept, r_value, p_value, std_err = stats.linregress(expShifts,
                                                                   calcShifts)
    scaled = [(x - intercept) / slope for x in calcShifts]

    return scaled


def BoltzmannWeight_DP5(Isomers,dp5Data):

    print("Conf", np.shape(np.array(dp5Data.ScaledAtomProbs)))

    for iso,Cscaled_probs,Cunscaled_probs,Hscaled_probs,Hunscaled_probs in zip( Isomers, dp5Data.CScaledAtomProbs,dp5Data.CUnscaledAtomProbs, dp5Data.HScaledAtomProbs,dp5DataHUnscaledAtomProbs):

        CB_scaled_probs = [0] * len(Cscaled_probs[0])

        CB_unscaled_probs = [0] * len(Cscaled_probs[0])

        HB_scaled_probs = [0] * len(Hscaled_probs[0])

        HB_unscaled_probs = [0] * len(Hscaled_probs[0])

        for population, Cconf_scaled_p,Cconf_unscaled_p, Hconf_scaled_p,Hconf_unscaled_p in zip(iso.Populations, Cscaled_probs,Cunscaled_probs, Hscaled_probs,Hunscaled_probs ):

            for i in range(len(CB_scaled_probs)):

                CB_scaled_probs[i] += Cconf_scaled_p[i] * population

                CB_unscaled_probs[i] += Cconf_unscaled_p[i] * population

            for i in range(len(HB_scaled_probs)):

                HB_scaled_probs[i] += Hconf_scaled_p[i] * population

                HB_unscaled_probs[i] += Hconf_unscaled_p[i] * population

        dp5Data.CBScaledAtomProbs.append(CB_scaled_probs)

        dp5Data.CBUnscaledAtomProbs.append(CB_unscaled_probs)

        dp5Data.HBScaledAtomProbs.append(HB_scaled_probs)

        dp5Data.HBUnscaledAtomProbs.append(HB_unscaled_probs)


    return dp5Data


def Calculate_DP5(dp5Data):

    for Cscaled_probs,Cunscaled_probs,Hscaled_probs,Hunscaled_probs in zip(dp5Data.CBScaledAtomProbs,dp5Data.CBUnscaledAtomProbs,dp5Data.HBScaledAtomProbs,dp5Data.CHBUnscaledAtomProbs):

        ## carbon

        CDP5unscaled = 1 - gmean([1 - p_si for p_si in Cunscaled_probs])

        CDP5scaled = (1 - gmean([1 - p_si for p_si in Cscaled_probs]))

        CDP5plus = (1  - gmean([1 - p_si for p_si in Cscaled_probs] + [1 - p_si for p_si in Cunscaled_probs]))

        dp5Data.CUnscaledprobs.append(CDP5unscaled)

        dp5Data.CScaledprobs.append(CDP5scaled)

        dp5Data.Cplusprobs.append(CDP5plus)

        ## proton

        HDP5unscaled = 1 - gmean([1 - p_si for p_si in Hunscaled_probs])

        HDP5scaled = (1 - gmean([1 - p_si for p_si in Hscaled_probs]))

        HDP5plus = (1 - gmean([1 - p_si for p_si in Hscaled_probs] + [1 - p_si for p_si in Hunscaled_probs]))

        dp5Data.HUnscaledprobs.append(HDP5unscaled)

        dp5Data.HScaledprobs.append(HDP5scaled)

        dp5Data.Hplusprobs.append(HDP5plus)

    return dp5Data


def Rescale_DP5(dp5Data,Settings):

    incorrect_kde = pickle.load(open(Path(Settings.ScriptDir) / "Ci_w_kde_mean_s_0.025.p" ,"rb"))

    correct_kde = pickle.load(open(Path(Settings.ScriptDir) / "Cc_w_kde_mean_s_0.025.p" ,"rb"))

    i = 0

    for scaled,unscaled in zip(dp5Data.CBScaledAtomProbs,dp5Data.CBUnscaledAtomProbs):

        dp5Data.CBScaledAtomProbs[i] = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in scaled  ]
        dp5Data.CBUnscaledAtomProbs[i] = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in unscaled  ]

        i += 1

    dp5Data.CDP5unscaledprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in dp5Data.CUnscaledprobs  ]  # Final DP5S
    dp5Data.CDP5scaledprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in dp5Data.CScaledprobs  ]   # Final DP5S
    dp5Data.CDP5plusprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in dp5Data.Cplusprobs  ]   # Final DP5s plus probs

    ###

    incorrect_kde = pickle.load(open(Path(Settings.ScriptDir) / "Hi_w_kde_mean_s_0.025.p" ,"rb"))

    correct_kde = pickle.load(open(Path(Settings.ScriptDir) / "Hc_w_kde_mean_s_0.025.p" ,"rb"))

    i = 0

    for scaled,unscaled in zip(dp5Data.HBScaledAtomProbs,dp5Data.HBUnscaledAtomProbs):

        dp5Data.HBScaledAtomProbs[i] = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in scaled  ]
        dp5Data.HBUnscaledAtomProbs[i] = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in unscaled  ]

        i += 1

    dp5Data.HDP5unscaledprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in dp5Data.HUnscaledprobs  ]  # Final DP5S
    dp5Data.HDP5scaledprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in dp5Data.HScaledprobs  ]   # Final DP5S
    dp5Data.HDP5plusprobs = [  (incorrect_kde.pdf(x) / (incorrect_kde.pdf(x) + correct_kde.pdf(x)))[0] for x in dp5Data.Hplusprobs  ]   # Final DP5s plus probs

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

                "Compounds": dp5Data.Compounds,
                "ConfCshifts": dp5Data.ConfCshifts,
                "AtomReps": dp5Data.AtomReps,
                "ConfHshifts": dp5Data.ConfCshifts,
                "HAtomReps": dp5Data.AtomReps,


                "CUnscaledAtomProbs": dp5Data.UnscaledAtomProbs,
                "CScaledAtomProbs": dp5Data.ScaledAtomProbs,
                "HUnscaledAtomProbs": dp5Data.UnscaledAtomProbs,
                "HScaledAtomProbs": dp5Data.ScaledAtomProbs,

                "CBUnscaledAtomProbs": dp5Data.BUnscaledAtomProbs,
                "CBScaledAtomProbs": dp5Data.BScaledAtomProbs,
                "HBUnscaledAtomProbs": dp5Data.BUnscaledAtomProbs,
                "HBScaledAtomProbs": dp5Data.BScaledAtomProbs,

                "CUnscaledprobs": dp5Data.CUnscaledprobs,
                "CScaledprobs": dp5Data.CScaledprobs,
                "Cplusprobs": dp5Data.Cplusprobs,
                "HUnscaledprobs": dp5Data.CUnscaledprobs,
                "HScaledprobs": dp5Data.CScaledprobs,
                "Cplusprobs": dp5Data.Cplusprobs,

                "CDP5unscaledprobs": dp5Data.DP5unscaledprobs,
                "CDP5scaledprobs": dp5Data.DP5scaledprobs,
                "CDP5plusprobs": dp5Data.DP5plusprobs,
                "HDP5unscaledprobs": dp5Data.DP5unscaledprobs,
                "HDP5scaledprobs": dp5Data.DP5scaledprobs,
                "HDP5plusprobs": dp5Data.DP5plusprobs}

    pickle.dump(data_dic , open(Path(Settings.OutputFolder) / "wf" / "data_dic.p","wb"))

    return dp5Data


def UnPickle_res(dp5Data,Settings):

    data_dic =  pickle.load(open(Path(Settings.OutputFolder) / "wf" / "data_dic.p","rb"))

    dp5Data.Cshifts = data_dic["Cshifts"]
    dp5Data.Cexp = data_dic["Cexp"]
    dp5Data.Clabels = data_dic["Clabels"]
    dp5Data.Hshifts = data_dic["Hshifts"]
    dp5Data.Hexp = data_dic["Hexp"]
    dp5Data.Hlabels = data_dic["Hlabels"]
    dp5Data.Cscaled = data_dic["Cscaled"]
    dp5Data.Hscaled = data_dic["Hscaled"]

    dp5Data.ConfCshifts = data_dic["ConfCshifts"]
    dp5Data.ConfHshifts = data_dic["ConfHshifts"]

    dp5Data.Compounds = data_dic["Compounds"]

    dp5Data.CAtomReps = data_dic["CAtomReps"]
    dp5Data.HAtomReps = data_dic["HAtomReps"]

    dp5Data.CUnscaledAtomProbs = data_dic["CUnscaledAtomProbs"]
    dp5Data.CScaledAtomProbs = data_dic["CScaledAtomProbs"]
    dp5Data.HUnscaledAtomProbs = data_dic["HUnscaledAtomProbs"]
    dp5Data.HScaledAtomProbs = data_dic["HScaledAtomProbs"]

    dp5Data.CBUnscaledAtomProbs = data_dic["CBUnscaledAtomProbs"]
    dp5Data.CBScaledAtomProbs = data_dic["CBScaledAtomProbs"]
    dp5Data.HBUnscaledAtomProbs = data_dic["HBUnscaledAtomProbs"]
    dp5Data.HBScaledAtomProbs = data_dic["HBScaledAtomProbs"]

    dp5Data.CUnscaledprobs = data_dic["CUnscaledprobs"]
    dp5Data.CScaledprobs = data_dic["CScaledprobs"]
    dp5Data.Cplusprobs = data_dic["Cplusprobs"]
    dp5Data.HUnscaledprobs = data_dic["HUnscaledprobs"]
    dp5Data.HScaledprobs = data_dic["HScaledprobs"]
    dp5Data.Hplusprobs = data_dic["Hplusprobs"]

    dp5Data.CDP5unscaledprobs = data_dic["CDP5unscaledprobs"]
    dp5Data.CDP5scaledprobs = data_dic["CDP5scaledprobs"]
    dp5Data.CDP5plusprobs = data_dic["CDP5plusprobs"]
    dp5Data.HDP5unscaledprobs = data_dic["HDP5unscaledprobs"]
    dp5Data.HDP5scaledprobs = data_dic["HDP5scaledprobs"]
    dp5Data.HDP5plusprobs = data_dic["HDP5plusprobs"]

    return dp5Data


def PrintAssignment(dp5Data):
    isomer = 0

    for Clabels, Cshifts, Cexp, Cscaled, atom_p in zip(dp5Data.Clabels, dp5Data.Cshifts, dp5Data.Cexp, dp5Data.Cscaled,dp5Data.CBScaledAtomProbs):
        dp5Data.output += ("\n\nAssigned C shifts for isomer " + str(isomer + 1) + ": ")

        PrintNMR(Clabels, Cshifts, Cscaled, Cexp,atom_p, dp5Data)

        isomer += 1

    for Hlabels, Hshifts, Hexp, Hscaled, atom_p in zip(dp5Data.Hlabels, dp5Data.Hshifts, dp5Data.Hexp, dp5Data.Hscaled,dp5Data.HBScaledAtomProbs):
        dp5Data.output += ("\n\nAssigned H shifts for isomer " + str(isomer + 1) + ": ")

        PrintNMR(Hlabels, Hshifts, Hscaled, Hexp,atom_p, dp5Data)

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
                           format(atom_p[i], "6.2f"))


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

    #### carbon

    c = 1

    for i in Isomers:
        dp5Data.output += "\nNumber of conformers for isomer " + str(c) + " = " + str(len(i.Conformers))

        c += 1

    PrintAssignment(dp5Data)

    dp5Data.output += ("\n\nResults of DP4 using Carbon: ")

    for i, p in enumerate(dp5Data.CDP5scaledprobs):
        dp5Data.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    dp5Data.output += ("\n\nResults of DP4: ")


    print("number of c carbons = " + str(len(Isomers[0].Clabels)))
    print("number of e carbons = " + str(len(dp5Data.Cexp[0])))

    print(dp5Data.output)


    ####proton

    dp5Data.output += ("\n\nResults of DP4 using Proton: ")

    for i, p in enumerate(dp5Data.HDP5scaledprobs):
        dp5Data.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    dp5Data.output += ("\n\nResults of DP4: ")


    print("number of c protons = " + str(len(Isomers[0].Hlabels)))
    print("number of e protons = " + str(len(dp5Data.Hexp[0])))

    print(dp5Data.output)


    if Settings.OutputFolder == '':

        out = open(str(os.getcwd()) + "/" + str(Settings.InputFiles[0] + "NMR.wf"), "w+")

    else:

        out = open(os.path.join(Settings.OutputFolder, str(Settings.InputFiles[0] + "NMR.wf")), "w+")

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