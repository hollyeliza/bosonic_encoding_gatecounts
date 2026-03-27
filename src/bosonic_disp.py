import numpy as np

def bosonic_disp_operator_matrix(d: int) -> np.ndarray:
    """
    This makes the matrix representation of the bosonic displacement operator 
    q = (a† + a) / sqrt(2) in a truncated Fock basis. Each element corresponds to the 
    matrix coefficient (strength of the transition) between the basis. Note that although 
    we encode the basis (unary, gray, sb) so that each state is represented using binary 
    (more than just the phonon occupation number).
    
    Note: The 1/sqrt.2 is needed to satisfy the commutation relation between position and 
    momentum - the commutator needs to be equal to i (where we set hbar to 1). The
    Liu et al. paper does not have the square root though.
    
    Note: This is the operator for position (displacement) if we consider the oscillations
    to be perfectly harmonic. The operators were chosen so that they satisfy the 
    Harmonic oscillator Hamiltonian and ensure the correct commutation relation.
 
    Parameters
    ----------
    d : int
        The highest occupation of the bosonic mode we want the physical system access to

    Returns
    -------
    np.ndarray
        (d × d) matrix representation of the position operator.
        
        Matrix containing the transition amplitudes between Fock states. The off-diagonal 
        elements scale as sqrt(n+1), and so as the occupation number increases it becomes
        'easier' to transition to even higher occupation states.
    """

    q = np.zeros((d, d), dtype=float) # first initialize an empty matrix of the right dimension
    for n in range(d - 1):
        val = np.sqrt(n + 1) / np.sqrt(2.0) # incorporating the square root into each value that is for normalization or smth
        q[n, n + 1] = val # adding this in 2 positions - matrix is tridiagonal
        q[n + 1, n] = val
    return q