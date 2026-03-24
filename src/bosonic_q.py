import numpy as np


def bosonic_q_matrix(d: int) -> np.ndarray:
    """
    This makes the matrix representation of the bosonic position operator 
    q = (a^\dagger + a)/sqrt(2) in the Fock basis {|0>, ..., |d-1> (as in
    the occupation, non encoded basis). The square root has been kept although
    not in the Liu et al. paper to seem everything the same to reproduce the
    Sawaya plots
    """
    q = np.zeros((d, d), dtype=float) # first initialize an empty matrix of the right dimension
    for n in range(d - 1):
        val = np.sqrt(n + 1) / np.sqrt(2.0) # incorporating the square root into each value that is for normalization or smth
        q[n, n + 1] = val # adding this in 2 positions - matrix is tridiagonal
        q[n + 1, n] = val
    return q

# displacement_operator_6 = bosonic_q_matrix(6)
# print(displacement_operator_6)

