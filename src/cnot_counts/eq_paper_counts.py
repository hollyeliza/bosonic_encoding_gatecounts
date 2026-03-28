# This code calculates the number of CNOTs using the equation in the paper.
# This assumes each Pauli string is implemented independently using a
# standard CNOT ladder with cost 2(p - 1). It is an overestimate as it doesn't
# account for any optimization (merging / cancelling gates). I just can't
# recreate the plots in the paper and trying to figure out what was done.

# It is an overestimate

# For a Hermitian pair |l><l'| + |l'><l|, the number of Pauli strings of length p is:

#     f(p; d_H, K) = (1/2) * 2^{d_H} * binom(K - d_H, p - d_H)

# where:
#     K   = |C_enc(l) ∪ C_enc(l')|
#     d_H = Hamming distance between encoded bitstrings R(l) and R(l')

# The naive upper bound on the CNOT count is then:

#     n_CNOT = sum_{p=2}^{K} f(p; d_H, K) * (2p - 2)

#            = sum_{p=2}^{K} (1/2) * 2^{d_H} * binom(K - d_H, p - d_H) * (2p - 2)

import math
from math import comb
from pauli_string_formation.encodings_b import bits_for_level, bitmask_subset

def hamming_distance(bits1: list[int], bits2: list[int]) -> int:
    """
    Return the Hamming distance between two equal-length bitstrings.
    """
    if len(bits1) != len(bits2):
        raise ValueError("Bitstrings must have the same length.")
    return sum(b1 != b2 for b1, b2 in zip(bits1, bits2))


def relevant_K(l: int, lp: int, d: int, encoding: str) -> int:
    """
    Return K = |C_enc(l) ∪ C_enc(lp)| from the Sawaya SI.
    """
    support_l = bitmask_subset(l, d, encoding)
    support_lp = bitmask_subset(lp, d, encoding)
    return len(support_l | support_lp)


def relevant_dH(l: int, lp: int, d: int, encoding: str) -> int:
    """
    Return d_H(R(l), R(lp)), the Hamming distance between the encoded levels.
    """
    bits_l = bits_for_level(l, d, encoding)
    bits_lp = bits_for_level(lp, d, encoding)
    return hamming_distance(bits_l, bits_lp)


def f_pauli_length_count(p: int, dH: int, K: int) -> int | float:
    """
    Sawaya SI Eq. (3): number of length-p Pauli strings for one Hermitian pair
    alpha|l><lp| + alpha*|lp><l|.

    f(p; dH, K) = (1/2) * 2^dH * C(K-dH, p-dH)

    Returns 0 if p is outside the allowed range.
    """
    if p < dH or p > K:
        return 0
    return 0.5 * (2 ** dH) * comb(K - dH, p - dH)

# comb (combinations): how many ways you can choose k items from n, without caring about order


def naive_cnot_upper_bound_from_dH_K(dH: int, K: int) -> int:
    """
    Sawaya SI Eq. (4): naive upper bound on CNOT count for one Hermitian pair
    alpha|l><lp| + alpha*|lp><l|, using 2(p-1) CNOTs for each length-p Pauli string.
    """
    total = 0
    for p in range(2, K + 1):
        total += f_pauli_length_count(p, dH, K) * (2 * p - 2)
    return int(round(total))


def naive_cnot_upper_bound_closed_form(dH: int, K: int) -> int:
    """
    Closed forms from Sawaya SI Eq. (5) for dH = 0, 1, 2.
    """
    if dH == 0:
        return (2 ** K) * K - 2 * (2 ** K) - 2
    elif dH == 1:
        return ((2 ** K) * K - (2 ** K)) // 2
    elif dH == 2:
        return ((2 ** K) * K) // 2
    else:
        raise ValueError("Closed form only given in the paper for dH = 0, 1, 2.")


def naive_cnot_upper_bound_for_pair(l: int, lp: int, d: int, encoding: str) -> int:
    """
    Compute the Sawaya naive upper-bound CNOT count for the Hermitian pair
    |l><lp| + |lp><l| in the chosen encoding.

    Uses:
      - K = |C_enc(l) ∪ C_enc(lp)|
      - dH = Hamming distance between encoded bitstrings
      - Eq. (4) from the SI
    """
    K = relevant_K(l, lp, d, encoding)
    dH = relevant_dH(l, lp, d, encoding)
    return naive_cnot_upper_bound_from_dH_K(dH, K)


def diagonal_naive_cnot_upper_bound(l: int, d: int, encoding: str) -> int:
    """
    Naive upper-bound CNOT count for a diagonal term |l><l|.

    This corresponds to dH = 0 and K = |C_enc(l)|.
    """
    K = len(bitmask_subset(l, d, encoding))
    return naive_cnot_upper_bound_from_dH_K(dH=0, K=K)