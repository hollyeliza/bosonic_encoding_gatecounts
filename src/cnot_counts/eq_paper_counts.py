# This code calculates the number of CNOTs using the equation in the paper.
# This assumes each Pauli string is implemented independently using a
# standard CNOT ladder with cost 2(p - 1). It is an overestimate as it doesn't
# account for any optimization (merging / cancelling gates).

# For a Hermitian pair |l><l'| + |l'><l|, the number of Pauli strings of length p is:

#     f(p; d_H, K) = (1/2) * 2^{d_H} * binom(K - d_H, p - d_H)

# where:
#     K   = |C_enc(l) ∪ C_enc(l')|
#     d_H = Hamming distance between encoded bitstrings R(l) and R(l')

# The naive upper bound on the CNOT count is then:

#     n_CNOT = sum_{p=2}^{K} f(p; d_H, K) * (2p - 2)

#            = sum_{p=2}^{K} (1/2) * 2^{d_H} * binom(K - d_H, p - d_H) * (2(p - 1))

from math import comb
from src.pauli_string_formation.encodings_b import bits_for_level, bitmask_subset


def hamming_distance(bits1: list[int], bits2: list[int]) -> int:
    """
    Return the Hamming distance between two equal-length bitstrings.

    The Hamming distance is the number of positions where the two bitstrings
    differ. This is used to compare the encoded forms R(l) and R(lp).
    """
    if len(bits1) != len(bits2):
        raise ValueError("Bitstrings must have the same length.")

    distance = 0
    for b1, b2 in zip(bits1, bits2):
        if b1 != b2:
            distance += 1

    return distance


def relevant_K(l: int, lp: int, d: int, encoding: str) -> int:
    """
    Return K = |C_enc(l) union C_enc(lp)|.

    K is the number of qubits in the union of the relevant bitmask subsets
    for the two encoded levels. For compact encodings such as standard binary
    and Gray, this is usually the full number of encoding qubits. For unary,
    K is smaller: usually 1 for diagonal terms and 2 for off-diagonal terms.
    """
    support_l = bitmask_subset(l, d, encoding)
    support_lp = bitmask_subset(lp, d, encoding)

    return len(support_l | support_lp)


def relevant_dH(l: int, lp: int, d: int, encoding: str) -> int:
    """
    Return the Hamming distance dH between the encoded bitstrings R(l) and R(lp).

    The Hamming distance counts the number of qubit positions where the two
    encoded levels differ. The bitstrings are represented in little-endian
    order, but since this convention is used consistently for both states,
    it does not affect the result.
    """
    bits_l = bits_for_level(l, d, encoding)
    bits_lp = bits_for_level(lp, d, encoding)

    return hamming_distance(bits_l, bits_lp)

# print(f'Hamming distance for |13><14| in Gray (expect 2):{relevant_dH(13, 14, 16, "sb")}')

def num_pauli_terms_length_p_diagonal(p: int, dH: int, K: int) -> int:
    """
    Return the number of length-p Pauli strings for a diagonal matrix element.

    For diagonal terms, no Hermitian conjugate is added. In this case dH should
    be 0, but the function keeps dH as an argument for consistency with the
    Sawaya notation.

    Returns 0 if p is outside the allowed range.

    Note: Un-needed for bosonic displacement (for completeness)
    """
    if p < dH or p > K:
        return 0

    return comb(K , p)


def num_pauli_terms_length_p_off_diagonal(p: int, dH: int, K: int) -> float:
    """
    Return the number of length-p Pauli strings for one off-diagonal Hermitian pair.

    This corresponds to a pair of transitions |l><lp| and |lp><l|. The count
    depends on the Hamming distance dH between the encoded levels and on K,
    the number of relevant qubits in the bitmask union.

    Returns 0 if p is outside the allowed range.
    """
    if p < dH or p > K:
        return 0.0

    return (2 ** dH) * comb(K - dH, p - dH)


def off_diagonal_naive_cnot_upper_bound_from_dH_K(dH: int, K: int) -> int:
    """
    Return the naive CNOT upper bound for one off-diagonal Hermitian pair.

    Each length-p Pauli string is assumed to be implemented independently with
    a standard CNOT ladder costing 2(p - 1) CNOTs. This does not include circuit
    optimization, cancellation, or merging between adjacent Pauli strings.

    The 0.5 from the Hermitian conjugate is counted here,
    """
    total = 0.0

    for p in range(2, K + 1):
        n_terms = num_pauli_terms_length_p_off_diagonal(p, dH, K)
        total += n_terms * (2 * (p - 1))

    return int(round(0.5 * total))


def diagonal_naive_cnot_upper_bound_from_dH_K(dH: int, K: int) -> int:
    """
    Return the naive CNOT upper bound for one diagonal matrix element.

    Each length-p Pauli string is assumed to be implemented independently with
    a standard CNOT ladder costing 2(p - 1) CNOTs. This is included for
    completeness, although the bosonic displacement operator has zero diagonal.

    Note: Un-needed for bosonic displacement (for completeness)
    """
    total = 0

    for p in range(2, K + 1):
        n_terms = num_pauli_terms_length_p_diagonal(p, dH, K)
        total += n_terms * (2 * (p - 1))

    return int(round(total))


def naive_cnot_upper_bound_from_dH_K(dH: int, K: int, diagonal: bool = False) -> int:
    """
    Return the naive CNOT upper bound for one matrix element contribution.

    If diagonal is True, the term is treated as a single diagonal matrix element.
    If diagonal is False, the term is treated as an off-diagonal Hermitian pair.
    """
    if diagonal: # if this is true l == lp
        return diagonal_naive_cnot_upper_bound_from_dH_K(dH, K)

    return off_diagonal_naive_cnot_upper_bound_from_dH_K(dH, K)


def naive_cnot_upper_bound_for_transition(
    l: int, lp: int, d: int, encoding: str
) -> int:
    """
    Return the naive CNOT upper bound for one transition |l><lp|.

    The function computes K and dH from the chosen encoding. Diagonal terms are
    counted once, while off-diagonal terms are treated as Hermitian pairs.
    """
    K = relevant_K(l, lp, d, encoding)
    dH = relevant_dH(l, lp, d, encoding)

    return naive_cnot_upper_bound_from_dH_K(
        dH=dH,
        K=K,
        diagonal=(l == lp),
    )
