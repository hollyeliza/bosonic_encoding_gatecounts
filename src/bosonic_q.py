import numpy as np


def bosonic_q_matrix(d: int) -> np.ndarray:
    """
    Truncated bosonic position operator q = (a^\dagger + a)/sqrt(2)
    in the Fock basis {|0>, ..., |d-1>}.
    """
    q = np.zeros((d, d), dtype=float)
    for n in range(d - 1):
        val = np.sqrt(n + 1) / np.sqrt(2.0)
        q[n, n + 1] = val
        q[n + 1, n] = val
    return q