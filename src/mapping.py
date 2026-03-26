# We encode the bosonic state (phonon occupation) from |l> to |R(l)> where
# R(l) represents the encoding (gray, sb, unary). To apply the bosonic
# operators like the displacement operator, we need |R(l)><R(l')|.
# To get this we break this outer product into a tensor product of one
# qubit outer products, that each have a nice mpping to a Pauli gate
# As the qubit can be 0 or 1, there are only 4 possibilities
#
# Example:
#
# |5><6|
# Gray: |0111><0101| Note: sequential occupations differ by 1 qubit - local!
# |0111><0101| = |0><0| x|1><1| x |1><0| x |1><1|
# 
# Each of these one qubit outer products can be replaced in in terms of X,Y,Z
# operators and this forms a sum of Pauli strings


from itertools import product
from openfermion import QubitOperator
from encodings_b import bits_for_level, bitmask_subset, n_qubits

def one_qubit_map(x: int, xp: int):
    """
    Return the Pauli decomposition of the single-qubit operator |x><xp|.

    The result is a list of (pauli_label, coefficient), where
    pauli_label ∈ {"I", "X", "Y", "Z"}.

    This is the building block used to construct multi-qubit Pauli
    strings for encoded transitions |R(l)><R(l')|.
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

# print(one_qubit_map(1,0))


def matrix_element_to_qubit_operator(
    l: int, lp: int, coeff: complex, d: int, encoding: str
) -> QubitOperator:
    """
    Convert a single logical bosonic matrix element, coeff * |l><lp|, into
    an OpenFermion QubitOperator written as a sum of Pauli strings.

    What this function starts from
    ------------------------------
    In the truncated bosonic Fock basis {|0>, |1>, ..., |d-1>}, a bosonic
    operator can be written as a matrix:
        A = sum_{l,lp} A[l,lp] |l><lp|.

    This function handles one such matrix element at a time:
        coeff * |l><lp|.

    The goal is to rewrite this logical transition as an operator acting on
    qubits, using the chosen encoding ("sb", "gray", or "unary").

    High-level idea
    ---------------
    1. Encode the two logical levels |l> and |lp> as qubit bitstrings:
           |l>  -> |R(l)>
           |lp> -> |R(lp)>.
    2. Rewrite the encoded outer product
           |R(l)><R(lp)|
       as a tensor product of single-qubit outer products:
           ⊗_q |x_q><x'_q|,
       where x_q and x'_q are the bits on qubit q.
    3. Replace each single-qubit outer product with its Pauli decomposition:
           |0><1| = 1/2 (X + iY) This corresponds to a single matrix that is a blend of these both
           |1><0| = 1/2 (X - iY)
           |0><0| = 1/2 (I + Z)
           |1><1| = 1/2 (I - Z).
    4. Multiply these local decompositions together to obtain a sum of
       multi-qubit Pauli strings.
    5. Return the result as a QubitOperator.

    Why the support is restricted
    -----------------------------
    For noncompact encodings such as unary, the encoded state does not occupy
    every qubit. In that case, the transition only needs to act on the qubits
    in C(l) ∪ C(lp), where C(l) is the qubit support of level l.

    For example, in unary:
        |5> -> |000001000...>
        |6> -> |000000100...>
    so the transition |5><6| only involves qubits 5 and 6.

    What the output is
    ------------------
    The returned object is an OpenFermion QubitOperator. This is a symbolic
    representation of a qubit operator as a sum of Pauli strings.

    For example, an operator such as
        0.25 * X5 X6 + 0.25j * X5 Y6 - 0.25j * Y5 X6 + 0.25 * Y5 Y6
    is stored as a QubitOperator with four separate Pauli terms.

    Initial misunderstanding: we are building a sum of multi-qubit operators not
    a sequence of gates 

    Concrete example
    ----------------
    Suppose:
        l = 5, lp = 6, coeff = 11, d = 12, encoding = "unary".

    Then
        |5><6|
    becomes a transition only on qubits 5 and 6:
        qubit 5: |1><0| = 1/2 (X - iY)
        qubit 6: |0><1| = 1/2 (X + iY).

    Therefore
        11 * |5><6|
    maps to
        11 * [1/2 (X - iY)]_5 ⊗ [1/2 (X + iY)]_6

    which expands to
        11/4 * ( X5 X6 + i X5 Y6 - i Y5 X6 + Y5 Y6 ).

    That full sum is what this function returns.

    Parameters
    ----------
    l : int
        Initial logical bosonic level.
    lp : int
        Final logical bosonic level.
    coeff : complex
        Coefficient multiplying |l><lp|.
    d : int
        Bosonic cutoff dimension.
    encoding : str
        Encoding used to map bosonic levels to qubits.
        Supported values: "sb", "gray", "unary".

    Returns
    -------
    QubitOperator
        The Pauli-string expansion of coeff * |l><lp| in the chosen encoding.
    """
    bits_l = bits_for_level(l, d, encoding)
    bits_lp = bits_for_level(lp, d, encoding)
    support = sorted(bitmask_subset(l, d, encoding) | bitmask_subset(lp, d, encoding)) # You do not need to act on every qubit if the encoding is noncompact.

    local_choices = []
    for q in support:
        local_choices.append(one_qubit_map(bits_l[q], bits_lp[q])) # a list of a list of pauli terms for each relevant qubit (in support)

    op = QubitOperator() # building operator
        # the operation on each qubit for a single transition is a single matrix (mix of X, Y, Z)
        # the tensor product is taken which results in Pauli strings (actual product of X, Y, Z correspondong to sequential application) 
        # This expands the tensor product. There are lots terms to mulitply out / combinations
        # I was caught up with the 'combination' but it is just how the X, Y, Z combine and need to, to multiply out tensor product
    for choice_tuple in product(*local_choices):
        term = []
        total_coeff = coeff
        for qubit, (pauli, c) in zip(support, choice_tuple):
            total_coeff *= c
            if pauli != "I":
                term.append((qubit, pauli))
        op += QubitOperator(tuple(term), total_coeff)
    return op

# print(f'In unary: {matrix_element_to_qubit_operator(5, 6, 7, 12, "unary")}') # only qubits 5 and 6 involved - more local as less in support 

# print(f'In gray: {matrix_element_to_qubit_operator(5, 6, 7, 12, "gray")}')


def matrix_to_qubit_operator(mat, d: int, encoding: str, tol: float = 1e-14) -> QubitOperator:
    """
    Takes a full dxd operator (logical matrix) and rewrites this as
    a sum of Pauli strings on each qubit.
    """
    op = QubitOperator()
    for l in range(d):
        for lp in range(d):
            coeff = mat[l, lp]
            if abs(coeff) > tol: # If the matrix element is numerically tiny, do nothing.
                op += matrix_element_to_qubit_operator(l, lp, coeff, d, encoding) # Convert that one matrix element into Pauli strings
    op.compress(abs_tol=1e-12) # combines duplicate strings and removes very tiny coefficients
    return op


def pauli_length(term) -> int:
    """"
    Returns the Pauli length of eaxh term in the string, 
    which is how many qubits this Pauli string acts on nontrivially.
    """
    return len(term)


def naive_cnot_count_from_qubit_operator(op: QubitOperator) -> int:
    """
    Count naive CNOTs term-by-term using 2(p-1) for each Pauli string of length p>1.
    Assumes each Pauli string is implemented separately with a standard CNOT ladder.
    p - 1 CNOTs to collect parity onto one qubit before applying Rz rotation and 
    reverse chain of another p-1 CNOTs -> total is 2(p-1) CNOTs per string.
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
    The amount of gate cancellation depends on the order of the terms...
    """
    pauli_priority = {"X": 0, "Y": 1, "Z": 2}

    def key_fn(item):
        term, coeff = item
        if len(term) == 0:
            return ((10_000, 10_000),) # practical trick making identity sort last with arbitrarily large number
        return tuple((q, pauli_priority[p]) for q, p in term)
    return sorted(op.terms.items(), key=key_fn)

# Note: QubitOperator is just a class