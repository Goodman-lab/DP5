import numpy as np
from numpy import linalg

'''
mol: array of coordinates
'''

def RMSDPrune(Isomers:list, settings):
    '''
    Takes a list of Isomer objects and removes redundant conformers

    Arguments:
    - Isomers: a list of Isomer objects. Must have Isomer.Conformers specified.
    - settings: Settings from PyDP4.py. Must have RMSD cutoff value (in Angstrom) and conformer limit specified

    Returns:
    - Isomers: a list of Isomer objects with modified list of conformers and the RMSD cutoff value
    '''
    for iso in Isomers:
        if len(iso.Conformers) < settings.PerStructConfLimit:
            continue
        else:
            conformers, cutoff = _StrictRMSDPrune(iso.Conformers, settings.InitialRMSDcutoff,
                                                 settings.PerStructConfLimit)
            iso.Conformers = conformers
            iso.RMSDCutoff = cutoff

    return Isomers

def _StrictRMSDPrune(conformers, cutoff, ConfLimit):
    '''
    Removes redundant conformers to fall below the specified conformer limit. Assumes conformers are sorted in ascending order by energy.

    Arguments:
    - conformers: a list of lists of lists of atomic coordinates. Follows the atom numbering from Isomer.InputFile.
    - atoms: a list of atomic symbols for atoms in Isomer. Follows the atom numbering from Isomer.InputFile. Redundant.
    - cutoff: a float. Conformers similar to lower energy conformers (RMSD<cutoff) are marked for deletion.
    - ConfLimit: an integer. Number of conformers should not exceed this limit.

    Returns:
    - PrunedConformers: a list of lists of lists of atomic coordinates. Follows the atom numbering from Isomer.InputFile.
    - cutoff: a float. The final value of the cutoff parameter.

    '''
    def confs_to_keep(rmsd_matrix:np.array, cutoff:float):
        '''
        Auxiliary function. Marks the non-redundant conformers. Assumes conformers are sorted in ascending order by energy.

        Arguments:
        - rmsd_matrix: an np.array containing aligned RMSD between conformers.
        - cutoff: a float. Conformers similar to lower energy conformers (RMSD<cutoff) are marked for deletion.
        Returns:
        - above_cutoff: a boolean array with indices set to True for retained conformers
        '''
        # checks for RMSDs above the cutoff
        mask = rmsd_matrix >= cutoff

        '''Sets the main diagonal and upper triangular elements as True. 
        Guarantees the retention of the redundant conformer encountered first.
        Higher energy redundant conformers should therefore have at least one False element in their column.'''
        mask[np.triu_indices_from(mask)]=True
        above_cutoff = mask.all(axis=1)

        return above_cutoff


    conf_arrays = [np.array(conformer,dtype=np.float32) for conformer in conformers]
    num_confs = len(conformers)
    #align them all

    conf_zeroed = []
    for conf in conf_arrays:
        conf_zeroed.append(_Move2Origin(conf))       

    rmsd_matrix = np.zeros((num_confs,num_confs))
    # populate the matrix
    for confid1 in range(0, num_confs):
        for confid2 in range(confid1+1,num_confs):
            res = _AlignedRMS(conf_zeroed[confid1], conf_zeroed[confid2])
            rmsd_matrix[confid1,confid2] = res
            rmsd_matrix[confid2,confid1] = res

    to_keep = confs_to_keep(rmsd_matrix, cutoff)

    while (to_keep.sum()) > ConfLimit:
        cutoff += 0.2
        to_keep = confs_to_keep(rmsd_matrix, cutoff)
    
    PrunedConformers = []
    for conf, keep in zip(conformers, to_keep):
        if keep:
            PrunedConformers.append(conf)

    return PrunedConformers, cutoff

def _AlignedRMS(mol1,mol2):
    '''
    aligns two conformers and estimates their RMSD. Does not take symmetry into account. Assumes conformers are translated to origin.

    arguments:
    - mol1: np.array of shape (3, number of atoms) containing coordinates for the first conformer
    - mol2: np.array of shape (3, number of atoms) containing coordinates for the second conformer
    - w: weighting vector. Previously used for weighted RMSD, now all are weighted equally.

    returns:
        - root mean squared deviation
    '''

    # align

    u = _qtrfit(mol1, mol2)
    mol1 = np.matmul(u, mol1.T).T

    # calculate RMS

    return _RMSMol(mol1,mol2)

def _RMSMol(mol1, mol2):
    '''
    computes RMS deviation of conformers. Assumes same ordering of the atoms.

    arguments:
    - mol1: np.array of shape (3, number of atoms) containing coordinates for the first conformer
    - mol2: np.array of shape (3, number of atoms) containing coordinates for the second conformer
    - w: weighting vector. Previously used for weighted RMSD, now all are weighted equally.

    returns:
        - root mean squared deviation
    '''
    # compute squared error
    sd = (mol1-mol2)**2

    return np.sqrt(sd.sum()/mol1.shape[0])

"""
 * This routine is a python translation of the QTRFIT Fortran routine
 * of David Heisterberg
 *
 * David J. Heisterberg
 * The Ohio Supercomputer Center, 1224 Kinnear Rd,
 * Columbus, OH  43212-1163
 * (614)292-6036; djh@osc.edu
 *
 * Python translation by K Ermanis 2015
 *
 * QTRFIT
 * Find the quaternion, q,[and left rotation matrix, u] 
 * that minimizes
 *   |qTBq - A| ^ 2   [|uB - A| ^ 2]
 * This is equivalent to maximizing Re (qTBTqA).
 * This is equivalent to finding the largest 
 * eigenvalue and corresponding
 * eigenvector of the matrix
 *  [V2   VUx  VUy  VUz ]
 *  [VUx  Ux2  UxUy UzUx]
 *  [VUy  UxUy Uy2  UyUz]
 *  [VUz  UzUx UyUz Uz2 ]
 * where
 *   V2   = Bx Ax + By Ay + Bz Az
 *   Ux2  = Bx Ax - By Ay - Bz Az
 *   Uy2  = By Ay - Bz Az - Bx Ax
 *   Uz2  = Bz Az - Bx Ax - By Ay
 *   VUx  = Bz Ay - By Az
 *   VUy  = Bx Az - Bz Ax
 *   VUz  = By Ax - Bx Ay
 *   UxUy = Bx Ay + By Ax
 *   UyUz = By Az + Bz Ay
 *   UzUx = Bz Ax + Bx Az
 * The left rotation matrix, u, is obtained from q by
 *   u = qT1q
 * INPUT
 *   n      - number of points
 *   b      - molecule to be rotated
 *   a      - reference molecule
 *   w      - weight vector
 * OUTPUT
 *   z      - eigenvalue
 *   q      - the best-fit quaternion
 *   u      - the best-fit left rotation matrix
 *   nr     - number of jacobi sweeps required
 
 KE:
 structure of molecule data:
 a[atomnr][coordinate]
 """

def _qtrfit(mol, refmol):
    '''
    generates quaternion rotation matrix for mol such as its deviation from a reference molecule is minimised
    '''
    # computes outer products of coordinates via a numpy trick, a being mol length, b, and c being coordinates
    product = np.einsum('ab,ac -> abc',mol,refmol)
    # outputs values
    xxyx, xxyy, xxyz, xyyx, xyyy, xyyz, xzyx, xzyy,xzyz = product.sum(axis=0).flatten()

    c = np.zeros((4,4),dtype=np.float64)
    c[0,0] = xxyx + xyyy + xzyz
    c[0,1] = c[1,0] = xzyy - xyyz
    c[1,1] = xxyx - xyyy - xzyz
    c[0,2] = c[2,0] = xxyz - xzyx
    c[1,2] = c[2,1] =  xxyy + xyyx
    c[2,2] = xyyy - xzyz - xxyx
    c[0,3] = c[3,0] = xyyx - xxyy
    c[1,3] = c[3,1] = xzyx + xxyz
    c[2,3] = c[3,2] = xyyz + xzyy
    c[3,3] = xzyz - xxyx - xyyy

    
    eval, evec = linalg.eigh(c)

    u = _q2mat(evec[3])

    return u

def _Move2Origin(mol:np.array):

    shift = mol.sum(axis=0)/mol.shape[0]
    mol = mol-shift

    return mol


def _q2mat(q:np.array):
    u = np.zeros((3,3),dtype=np.float64)

    u[0,0] = q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2
    u[1,0] = 2.0 *(q[1] * q[2] - q[0] * q[3])
    u[2,0] = 2.0 *(q[1] * q[3] + q[0] * q[2])
    u[0,1] = 2.0 *(q[2] * q[1] + q[0] * q[3])
    u[1,1] = q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2
    u[2,1] = 2.0 *(q[2] * q[3] - q[0] * q[1])
    u[0,2] = 2.0 *(q[3] * q[1] - q[0] * q[2])
    u[1,2] = 2.0 *(q[3] * q[2] + q[0] * q[1])
    u[2,2] = q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2
    
    return u
