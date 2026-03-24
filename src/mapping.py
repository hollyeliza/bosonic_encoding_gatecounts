from itertools import product
from openfermion import QubitOperator
from encodings_b import bits_for_level, bitmask_subset, n_qubits


def one_qubit_map(x: int, xp: int):
    """
    Return local map for |x><xp| as a list of (pauli_label, coeff),
    where pauli_label in {"I","X","Y","Z"}.
    """
    if x == 0 and xp == 1:
        return [("X", 0.5), ("Y", 0.5j)]
    if x == 1 and xp == 0:
        return [("X", 0.5), ("Y", -0.5j)]
    if x == 0 and xp == 0:
        return [("I", 0.5), ("Z", 0.5)]
    if x == 1 and xp == 1:
        return [("I", 0.5), ("Z", -0.5)]
    raise ValueError("Bits must be 0 or 1.")


def matrix_element_to_qubit_operator(
    l: int, lp: int, coeff: complex, d: int, encoding: str
) -> QubitOperator:
    """
    Map coeff * |l><lp| to a QubitOperator using the Sawaya rules.
    For noncompact encodings, only use C(l) U C(lp).
    """
    bits_l = bits_for_level(l, d, encoding)
    bits_lp = bits_for_level(lp, d, encoding)
    support = sorted(bitmask_subset(l, d, encoding) | bitmask_subset(lp, d, encoding))

    local_choices = []
    for q in support:
        local_choices.append(one_qubit_map(bits_l[q], bits_lp[q]))

    op = QubitOperator()
    for choice_tuple in product(*local_choices):
        term = []
        total_coeff = coeff
        for qubit, (pauli, c) in zip(support, choice_tuple):
            total_coeff *= c
            if pauli != "I":
                term.append((qubit, pauli))
        op += QubitOperator(tuple(term), total_coeff)
    return op


def matrix_to_qubit_operator(mat, d: int, encoding: str, tol: float = 1e-14) -> QubitOperator:
    """
    Map a full dxd operator to qubits by summing all nonzero entries.
    """
    op = QubitOperator()
    for l in range(d):
        for lp in range(d):
            coeff = mat[l, lp]
            if abs(coeff) > tol:
                op += matrix_element_to_qubit_operator(l, lp, coeff, d, encoding)
    op.compress(abs_tol=1e-12)
    return op


def pauli_length(term) -> int:
    return len(term)


def naive_cnot_count_from_qubit_operator(op: QubitOperator) -> int:
    """
    Count naive CNOTs term-by-term using 2(p-1) for each Pauli string of length p>1.
    """
    total = 0
    for term, coeff in op.terms.items():
        if abs(coeff) < 1e-12:
            continue
        p = pauli_length(term)
        if p > 1:
            total += 2 * (p - 1)
    return total


def sorted_terms_pseudo_alphabetical(op: QubitOperator):
    """
    Try to mimic the paper's pseudo-alphabetical ordering:
    X on low qubits first, then Y, then Z, etc.
    """
    pauli_priority = {"X": 0, "Y": 1, "Z": 2}

    def key_fn(item):
        term, coeff = item
        if len(term) == 0:
            return ((10_000, 10_000),)
        return tuple((q, pauli_priority[p]) for q, p in term)

    return sorted(op.terms.items(), key=key_fn)