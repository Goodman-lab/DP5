# DP5


[DP5](https://doi.org/10.1039/D1SC04406K) integrates NMR-AI, software for automatic processing, assignment
and visualisation of raw NMR data. This functionality affords fully automated DP4 analysis of databases of molecules
with no user input. DP5 also utilises a fully revised version of PyDP4 updated for Python 3 and with major workflow
improvements.

## Latest release

This latest release represents a leap forward in structural uncertainty calculation through the addition of the new DP5
module. By adding this functionality to our existing state of the art DP5 software, both DP4 and DP5
probabilities can be calculated alongside each other fully automatically.

DP5 as described in our [latest paper](https://doi.org/10.1039/D1SC04406K) calculates the standalone probability of
structures being correct.  If a user cannot guarantee the correct structure is in the list of proposals, DP5 calculates
the probability each is correct, quantifying this uncertainty.

When the user is certain a one of their proposals is correct i.e in the case of elucidating relative stereochemistry,
DP4 is the best metric to use, as this uses the additional information that one of the proposals must be correct.

The power of DP5 and DP4 however is greatest when they are used in conjunction as DP5 can be used to increase
reliability and confidence in the conclusions of DP4 calculations (also described in the paper).


## Graphical user interface

Users can now also choose to run DP5 within a Graphical user interface for increased ease of use when single DP4
calculations. This graphical user interface can also be used to explore the spectra processed and assigned by NMR-AI.
DP5 can also be utilised to perform assignments of NMR spectra for a single isomer of a molecule. These assignments
can be viewed in the graphical output from DP5 or investigated interactively using The included GUI application.

## Further details

More details about the software can be found in the publication
[DP5 Straight from Spectrometer to Structure](https://doi.org/10.1039/D1SC04406K) and in the accompanying supporting
information.
