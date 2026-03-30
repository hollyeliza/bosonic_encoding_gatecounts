
from itertools import product
from openfermion import QubitOperator
from src.pauli_string_formation.encodings_b import bits_for_level, bitmask_subset, n_qubits
from src.pauli_string_formation.bosonic_disp import bosonic_disp_operator_matrix

def one_qubit_map(x: int, xp: int):
    """
    Return the Pauli decomposition of the single-qubit operators:
    |x><xp| (|0><0|,|1><1|, |1><0|, |1><1|)

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


def single_matrix_element_to_qubit_operator(
    l: int, lp: int, coeff: complex, d: int, encoding: str
) -> QubitOperator:
    """
    Builds a QubitOperator object for a given transition - in other words converts a 
    single bosonic matrix element coeff * |l><lp| into a sum of Pauli strings acting on 
    qubits (we need to do this for all transitions to get a big string of Pauli strings!).

    'single' as h.c. logic in next function.

    A bosonic operator written in the truncated Fock basis can be decomposed as
    a sum over matrix elements |l><lp|. This function takes one such element and
    maps it into qubit space using the chosen encoding (standard binary, Gray, or unary).

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
        
        From QubitOperator class docstring:  A QubitOperator represents a sum 
        of terms acting on qubits and overloads operations for easy manipulation 
        of these objects by the user.
    """
    bits_l = bits_for_level(l, d, encoding)
    bits_lp = bits_for_level(lp, d, encoding)
    # You do not need to act on every qubit if the encoding is noncompact:
    support = sorted(bitmask_subset(l, d, encoding) | bitmask_subset(lp, d, encoding)) 
    local_choices = [] # the unique, sorted qubit numbers of relvant qubits
    for q in support:
        local_choices.append(one_qubit_map(bits_l[q], bits_lp[q])) # a list of a list of pauli terms for each relevant qubit (in support) 
    # print(local_choices)
    # [[('X', 0.5), ('Y', (-0-0.5j))], [('X', 0.5), ('Y', 0.5j)]] This is the one qubit operation on 5 and 6 (rest identity)

    op = QubitOperator() # A sum of Pauli strings represented nicely
    # making an object of type QubitOperator to store the Pauli string which has many attributes
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

def matrix_element_to_qubit_operator(
    l: int, lp: int, coeff: complex, d: int, encoding: str
) -> QubitOperator:
    """
    Builds a QubitOperator object for a given transition - in other words converts a 
    single bosonic matrix element coeff * |l><lp| into a sum of Pauli strings acting on 
    qubits (we need to do this for all transitions to get a big string of Pauli strings!).

    A bosonic operator written in the truncated Fock basis can be decomposed as
    a sum over matrix elements |l><lp|. This function takes one such element and
    maps it into qubit space using the chosen encoding (standard binary, Gray, or unary).

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
        
        From QubitOperator class docstring:  A QubitOperator represents a sum 
        of terms acting on qubits and overloads operations for easy manipulation 
        of these objects by the user.
    """

    op_forward = single_matrix_element_to_qubit_operator(l, lp, coeff, d, encoding)
    op_backward = single_matrix_element_to_qubit_operator(lp, l, coeff.conjugate(), d, encoding)

    op = op_forward + op_backward
    op.compress(abs_tol=1e-12)
    return op


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
                # print(f'next: {op}')
                op += matrix_element_to_qubit_operator(l, lp, coeff, d, encoding) # Convert that one matrix element into Pauli strings
    op.compress(abs_tol=1e-12) # combines duplicate strings and removes very tiny coefficients
    return op


def pseudo_alphabetical_qubit_operator(op: QubitOperator) -> QubitOperator:
    """
    This takes in the QubitOperator Pauli string from the matrix_to_qubit_operator
    (in other words the bosonic displacement operator in Pauli strings) and 
    returns a new Pauli string whose terms are inserted in the "pseudo-alphabetical" order
    that was described in the paper.

    Returns another QubitOperator (just ordered)
    """

    pauli_priority = {"X": 0, "Y": 1, "Z": 2}

    def key_fn(item):
        term, coeff = item

        if len(term) == 0:
            return (10_000,)

        # ensure sorted by qubit index
        term = sorted(term, key=lambda x: x[0])

        first_qubit, first_pauli = term[0]
        primary = pauli_priority[first_pauli]

        qubit_pattern = tuple(q for q, _ in term)
        pauli_pattern = tuple(pauli_priority[p] for _, p in term)

        return (primary, qubit_pattern, pauli_pattern)

    # sort terms
    sorted_items = sorted(op.terms.items(), key=key_fn)

    # rebuild operator in that order
    new_op = QubitOperator()
    for term, coeff in sorted_items:
        new_op += QubitOperator(term, coeff)

    return new_op
