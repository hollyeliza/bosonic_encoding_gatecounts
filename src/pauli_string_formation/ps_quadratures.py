import numpy as np

# The position (q) and momentum (p) operators are known as the phase space quadrature operators.
# They are expresse in terms of the bosonic ladder operators and generate displacements in phase
# space. The position operator generates displacements along the momentum axis, and the momentum 
# operator generates displacement along the position axis. This code focuses calculates the CNOT
# counts for the position operator matrix only, but can easily be modified for the momentum 
# operator case. The momentum operator matrix is therefore also below, but not used further.

# Additional notes:

# 1) This code generates the matrix representation of the position (and momentum) operators in
# a truncated Fock basis. Each element corresponds to the transition amplitudes between Fock states.
# Tha matrix is tri-diagonal with zeros on the diagonal. As the non-zero elements scale as sqrt(n+1), 
# as the occupation number increases it becomes 'easier' to transition to even higher occupation states 
# (Note: very highly occupied states are less important as they are unstable and have short life times 
# (Ella/Drew)).

# 2) The 1/sqrt.2 is needed to satisfy the commutation relation between position and 
# momentum - the commutator needs to be equal to i (where we set hbar to 1). The
 # Liu et al. paper does not have the square root though.
    
# 2) These are only the operators if we consider the oscillations to be perfectly harmonic. The 
# operators were chosen so that they satisfy the Harmonic oscillator Hamiltonian and ensure the 
# correct commutation relation.


def position_operator_matrix(d: int) -> np.ndarray:
    """
    Matrix representation of the position operator

        p = (a† + a) / sqrt(2)
    
    Parameters
    ----------
    d : int
        The highest occupation of the bosonic mode we want the physical system access to - 1
        ( d = N_max + 1 )

    Returns
    -------
    np.ndarray
        (d × d) matrix representation of the position operator.
    """

    q = np.zeros((d, d), dtype=float) # first initialize an empty matrix of the right dimension
    for n in range(d - 1):
        val = np.sqrt(n + 1) / np.sqrt(2.0) # incorporating the square root into each valuemath
        q[n, n + 1] = val # adding this in 2 positions - matrix is tridiagonal
        q[n + 1, n] = val
    return q

# print(bosonic_disp_operator_matrix(9))

# [[0.         0.70710678 0.         0.         0.         0.
#   0.         0.         0.        ]
#  [0.70710678 0.         1.         0.         0.         0.
#   0.         0.         0.        ]
#  [0.         1.         0.         1.22474487 0.         0.
#   0.         0.         0.        ]
#  [0.         0.         1.22474487 0.         1.41421356 0.
#   0.         0.         0.        ]
#  [0.         0.         0.         1.41421356 0.         1.58113883
#   0.         0.         0.        ]
#  [0.         0.         0.         0.         1.58113883 0.
#   1.73205081 0.         0.        ]
#  [0.         0.         0.         0.         0.         1.73205081
#   0.         1.87082869 0.        ]
#  [0.         0.         0.         0.         0.         0.
#   1.87082869 0.         2.        ]
#  [0.         0.         0.         0.         0.         0.
#   0.         2.         0.        ]]



def momentum_operator_matrix(d: int) -> np.ndarray:
    """
    Matrix representation of the momentum operator

        p = i(a† - a) / sqrt(2)

    Parameters
    ----------
    d : int
        The highest occupation of the bosonic mode we want the physical system access to - 1
        ( d = N_max + 1 )
    
    Returns
    -------
    np.ndarray
        (d × d) matrix representation of the position operator.
    """
    p = np.zeros((d, d), dtype=complex)

    for n in range(d - 1):
        val = np.sqrt(n + 1) / np.sqrt(2.0)

        p[n, n + 1] = -1j * val
        p[n + 1, n] = 1j * val

    return p


# print(bosonic_momentum_operator_matrix(9))

# [[0.+0.j         0.-0.70710678j 0.+0.j         0.+0.j
#   0.+0.j         0.+0.j         0.+0.j         0.+0.j
#   0.+0.j        ]
#  [0.+0.70710678j 0.+0.j         0.-1.j         0.+0.j
#   0.+0.j         0.+0.j         0.+0.j         0.+0.j
#   0.+0.j        ]
#  [0.+0.j         0.+1.j         0.+0.j         0.-1.22474487j
#   0.+0.j         0.+0.j         0.+0.j         0.+0.j
#   0.+0.j        ]
#  [0.+0.j         0.+0.j         0.+1.22474487j 0.+0.j
#   0.-1.41421356j 0.+0.j         0.+0.j         0.+0.j
#   0.+0.j        ]
#  [0.+0.j         0.+0.j         0.+0.j         0.+1.41421356j
#   0.+0.j         0.-1.58113883j 0.+0.j         0.+0.j
#   0.+0.j        ]
#  [0.+0.j         0.+0.j         0.+0.j         0.+0.j
#   0.+1.58113883j 0.+0.j         0.-1.73205081j 0.+0.j
#   0.+0.j        ]
#  [0.+0.j         0.+0.j         0.+0.j         0.+0.j
#   0.+0.j         0.+1.73205081j 0.+0.j         0.-1.87082869j
#   0.+0.j        ]
#  [0.+0.j         0.+0.j         0.+0.j         0.+0.j
#   0.+0.j         0.+0.j         0.+1.87082869j 0.+0.j
#   0.-2.j        ]
#  [0.+0.j         0.+0.j         0.+0.j         0.+0.j
#   0.+0.j         0.+0.j         0.+0.j         0.+2.j
#   0.+0.j        ]]