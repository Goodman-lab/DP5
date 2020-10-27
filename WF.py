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


c_distance = 4.532297920317418


class WFdata:
    def __init__(self,ScriptPath):


        self.Cshifts = []  # Carbon shifts used in DP4 calculation
        self.Cexp = []  # Carbon experimental shifts used in DP4 calculation
        self.Clabels = []  # Carbon atom labels
        self.Hshifts = []  # Proton shifts used in DP4 calculation
        self.Hexp = []  # Proton experimental shifts used in DP4 calculation
        self.Hlabels = []  # Proton atom labels
        self.Cscaled = []  # Internally scaled carbon shifts
        self.Hscaled = []  # Internally scaled proton shifts

        self.Cscalederrors = []  # Scaled Carbon prediction errors
        self.Hscalederrors = []  # Scaled Proton prediction errors

        self.Cerrors = []  # Uncaled Carbon prediction errors
        self.Herrors = []  # Unscaled Proton prediction errors


        self.Compounds = [] #qml compound objects for each isomer
        self.Atomreps = [] # FCHL representations ordered by Clabels

        self.Cprobs = []  # Scaled carbon prediction error probabilities
        self.Hprobs = []  # Scaled proton prediction error probabilities

        self.WFprobs = []  # Final WFS
        self.WFplusprobs = []  # Final WFs plus probs


        self.output = str()  # final DP4 output


        self.folded_scaled_errors = pickle.load(open(ScriptPath /  "folded_scaled_errors.p", "rb"))

        self.atomic_reps = pickle.load(open(ScriptPath /  "atomic_reps_c.p", "rb"))

        self.mean_abs_error = np.mean(abs(self.folded_scaled_errors))





def worry_factor(wfData,sigma):

    #for each atom in the molecule calculate the atomic worry factor

    for iso in range(len(wfData)):

        errors =  wfData.Cerrors[iso]

        scaled_errors =  wfData.Cscalederrors[iso]

        reps =  wfData.Atomreps[iso]

        probs = []
        scaled_probs = []

        for e,s_e,r in zip(errors,scaled_errors,reps):

            # calculate similarites between this atom and those in the atomic representation test set

            K_sim = get_atomic_kernels(np.array([r]), wfData.atomic_reps, [sigma],
                                       cut_distance=c_distance)[0][0]

            K_sim = np.hstack((K_sim, K_sim))

            # calculate kde using K_sim as the weighting function

            kde_estimator = kde(wfData.folded_scaled_errors, weights=K_sim)

            e_diff = abs(e - wfData.mean_abs_error)

            p = kde_estimator.integrate_box_1d(wfData.mean_abs_error - e_diff, wfData.mean_abs_error + e_diff)

            probs.append(p)

            s_e_diff = abs(s_e - wfData.mean_abs_error)

            s_p = kde_estimator.integrate_box_1d(wfData.mean_abs_error - s_e_diff, wfData.mean_abs_error + s_e_diff)

            scaled_probs.append(s_p)

        inverse_mean_wf = 1 - gmean([1 - p_si for p_si in probs])

        inverse_mean_wf_plus = (1 - gmean([1 - p_si for p_si in scaled_probs])) * inverse_mean_wf

        wfData.WFprobs.append(inverse_mean_wf)
        wfData.WFplusprobs.append(inverse_mean_wf_plus)

    return wfData


def ProcessIsomers(WFdata, Isomers,Settings):
    # extract calculated and experimental shifts and add to WFdata instance

    # Carbon

    # make sure any shifts with missing peaks are removed from all isomers

    removedC = []

    removedH = []

    for iso in Isomers:

        WFdata.Cshifts.append([])
        WFdata.Cexp.append([])
        WFdata.Clabels.append([])

        for shift, exp, label in zip(iso.Cshifts, iso.Cexp, iso.Clabels):

            if exp != '':
                WFdata.Cshifts[-1].append(shift)
                WFdata.Cexp[-1].append(exp)
                WFdata.Clabels[-1].append(label)
                WFdata.Cerrors[-1].append(abs(shift - exp))

            elif label not in removedC:

                removedC.append(label)

    for l in removedC:

        for j, Clabel in enumerate(WFdata.Clabels):

            if l in Clabel:
                i = Clabel.index(l)

                WFdata.Cshifts[j].pop(i)

                WFdata.Cexp[j].pop(i)

                WFdata.Clabels[j].pop(i)

    # proton
    for iso in Isomers:

        WFdata.Hshifts.append([])
        WFdata.Hexp.append([])
        WFdata.Hlabels.append([])

        for shift, exp, label in zip(iso.Hshifts, iso.Hexp, iso.Hlabels):

            if exp != '':
                WFdata.Hshifts[-1].append(shift)
                WFdata.Hexp[-1].append(exp)
                WFdata.Hlabels[-1].append(label)
                WFdata.Herrors[-1].append(abs(shift - exp))

            elif label not in removedH:

                removedH.append(label)

    for l in removedH:

        for j, Hlabel in enumerate(WFdata.Hlabels):

            if l in Hlabel:
                i = Hlabel.index(l)

                WFdata.Hshifts[j].pop(i)

                WFdata.Hexp[j].pop(i)

                WFdata.Hlabels[j].pop(i)

    #write qml compound objects and atomic representations

    #first have to write xyz file for each input structure

    for iso in Isomers:

        #open new xyz file

        xyz_file = open(Settings.OutputFolder / iso.InputFile.stem +  ".xyz")

        #find conformer with the lowest energy

        ind = iso.Populations.index(max(iso.Populations))

        geom = iso.DFTConformers[ind]

        xyz_file.write(str(len(iso.Atoms)) + "\n" + "\n")

        for atom, coords in zip(iso.Atoms, geom):

            xyz_file.write(atom + " " + str(coords[0]) + " " + str(coords[1]) + " " + str(coords[2]) + "\n")

        xyz_file.close()

        WFdata.Compounds.append(qml.Compound(xyz = Settings.OutputFolder / iso.InputFile.stem +  ".xyz"))

        WFdata.Atomreps.append([])

        for C_l in iso.Clabels:

            ind = int(C_l.split("C")[1])

            WFdata.Atomreps[-1].append(WFdata.Compounds[-1][ind])

    return WFdata


def InternalScaling(WFdata):
    # perform internal scaling process

    # calculate prediction errors

    if len(WFdata.Cexp[0]) > 0:

        for Cshifts, Cexp in zip(WFdata.Cshifts, WFdata.Cexp):
            WFdata.Cscaled.append(ScaleNMR(Cshifts, Cexp))

        for Cscaled, Cexp in zip(WFdata.Cscaled, WFdata.Cexp):
            WFdata.Cscalederrors.append([Cscaled[i] - Cexp[i] for i in range(0, len(Cscaled))])

    if len(WFdata.Hexp[0]) > 0:

        for Hshifts, Hexp in zip(WFdata.Hshifts, WFdata.Hexp):
            WFdata.Hscaled.append(ScaleNMR(Hshifts, Hexp))

        for Hscaled, Hexp in zip(WFdata.Hscaled, WFdata.Hexp):
            WFdata.Hscalederrors.append([Hscaled[i] - Hexp[i] for i in range(0, len(Hscaled))])

    return WFdata


def ScaleNMR(calcShifts, expShifts):
    slope, intercept, r_value, p_value, std_err = stats.linregress(expShifts,
                                                                   calcShifts)
    scaled = [(x - intercept) / slope for x in calcShifts]

    return scaled



def PrintAssignment(WFdata):
    isomer = 0

    for Clabels, Cshifts, Cexp, Cscaled in zip(WFdata.Clabels, WFdata.Cshifts, WFdata.Cexp, WFdata.Cscaled):
        WFdata.output += ("\n\nAssigned C shifts for isomer " + str(isomer + 1) + ": ")

        PrintNMR(Clabels, Cshifts, Cscaled, Cexp, WFdata)

        isomer += 1

    isomer = 0

    for Hlabels, Hshifts, Hexp, Hscaled in zip(WFdata.Hlabels, WFdata.Hshifts, WFdata.Hexp, WFdata.Hscaled):
        WFdata.output += ("\n\nAssigned H shifts for isomer " + str(isomer + 1) + ": ")

        PrintNMR(Hlabels, Hshifts, Hscaled, Hexp, WFdata)

        isomer += 1


def PrintNMR(labels, values, scaled, exp, WFdata):
    s = np.argsort(values)

    svalues = np.array(values)[s]

    slabels = np.array(labels)[s]

    sscaled = np.array(scaled)[s]

    sexp = np.array(exp)[s]

    WFdata.output += ("\nlabel, calc, corrected, exp, error")

    for i in range(len(labels)):
        WFdata.output += ("\n" + format(slabels[i], "6s") + ' ' + format(svalues[i], "6.2f") + ' '
                           + format(sscaled[i], "6.2f") + ' ' + format(sexp[i], "6.2f") + ' ' +
                           format(sexp[i] - sscaled[i], "6.2f"))


def MakeOutput(WFdata, Isomers, Settings):
    # add some info about the calculation

    WFdata.output += Settings.InputFiles[0] + "\n"

    WFdata.output += "\n" + "Solvent = " + Settings.Solvent

    WFdata.output += "\n" + "Force Field = " + Settings.ForceField + "\n"

    if 'o' in Settings.Workflow:
        WFdata.output += "\n" + "DFT optimisation Functional = " + Settings.oFunctional
        WFdata.output += "\n" + "DFT optimisation Basis = " + Settings.oBasisSet

    if 'e' in Settings.Workflow:
        WFdata.output += "\n" + "DFT energy Functional = " + Settings.eFunctional
        WFdata.output += "\n" + "DFT energy Basis = " + Settings.eBasisSet

    if 'n' in Settings.Workflow:
        WFdata.output += "\n" + "DFT NMR Functional = " + Settings.nFunctional
        WFdata.output += "\n" + "DFT NMR Basis = " + Settings.nBasisSet

    if Settings.StatsParamFile != "none":
        WFdata.output += "\n\nStats model = " + Settings.StatsParamFile

    WFdata.output += "\n\nNumber of isomers = " + str(len(Isomers))

    c = 1

    for i in Isomers:
        WFdata.output += "\nNumber of conformers for isomer " + str(c) + " = " + str(len(i.Conformers))

        c += 1

    PrintAssignment(WFdata)

    WFdata.output += ("\n\nResults of DP4 using Proton: ")

    for i, p in enumerate(WFdata.HDP4probs):
        WFdata.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    WFdata.output += ("\n\nResults of DP4 using Carbon: ")

    for i, p in enumerate(WFdata.CDP4probs):
        WFdata.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    WFdata.output += ("\n\nResults of DP4: ")

    for i, p in enumerate(WFdata.DP4probs):
        WFdata.output += ("\nIsomer " + str(i + 1) + ": " + format(p * 100, "4.1f") + "%")

    print("number of c protons = " + str(len(Isomers[0].Hlabels)))
    print("number of c carbons = " + str(len(Isomers[0].Clabels)))

    print("number of e protons = " + str(len(WFdata.Hexp[0])))
    print("number of e carbons = " + str(len(WFdata.Cexp[0])))

    print(WFdata.output)

    if Settings.OutputFolder == '':

        out = open(str(os.getcwd()) + "/" + str(Settings.InputFiles[0] + "NMR.dp4"), "w+")

    else:

        out = open(os.path.join(Settings.OutputFolder, str(Settings.InputFiles[0] + "NMR.dp4")), "w+")

    out.write(WFdata.output)

    out.close()

    return WFdata


#open and read the input geometery NMR data etc

#generate fchl representation for the new molecule

#open required lists for kde process

#do the worry factor

