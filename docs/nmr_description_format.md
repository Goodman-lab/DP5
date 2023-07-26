# NMR description format

!!! note TODO

    Add brief description of what the file stores/goals of the format etc.

## Structure

Sections are seperated by empty lines.

1. The first section is assigned C shifts, can also be (any).
2. Second section is (un)assigned H shifts.
3. This section defines chemically equivalent atoms. Each line is a new set,
   all atoms in a line are treated as equivalent, their computed shifts averaged.
4. Final section, starting with a keyword OMIT defines atoms to be ignored.
   Atoms defined in this section do not need a corresponding shift in the NMR
   description


## Example file

An example NMRFILE is as follows:

```
59.58(C3),127.88(C11),127.52(C10),115.71(C9),157.42(C8),133.98(C23),118.22(C22),115.79(C21),158.00(C20),167.33(C1),59.40(C2),24.50(C31),36.36(C34),71.05(C37),142.14(C42),127.50(C41),114.64(C40),161.02(C39)

4.81(H5),7.18(H15),6.76(H14),7.22(H28),7.13(H27),3.09(H4),1.73(H32 or H33),1.83(H32 or H33),1.73(H36 or H35),1.73(H36 or H35),4.50(H38),7.32(H47),7.11(H46)

H15,H16
H14,H17
H28,H29
H27,H30
H47,H48
H46,H49
C10,C12
C9,C13
C22,C24
C21,C25
C41,C43
C40,C44

OMIT H19,H51
```

