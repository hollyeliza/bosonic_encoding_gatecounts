from itertools import product
from openfermion import QubitOperator
from src.pauli_string_formation.encodings_b import bits_for_level, bitmask_subset

def one_qubit_map(x: int, xp: int):
    """
    Returns the Pauli decomposition of the single-qubit operators:
    |x><xp| (|0><0|,|1><1|, |1><0|, |1><1|)

    This is the building block used to construct multi-qubit Pauli
    strings for encoded transitions |R(l)><R(l')|.

    Parameters
    ----------
    x : int
        Computational-basis value of the ket state. Must be either 0 or 1.
    xp : int
        Computational-basis value of the bra state. Must be either 0 or 1.

    Returns
    -------
    list of tuples (pauli_label, coefficient), where 
    pauli_label ∈ {"I", "X", "Y", "Z"}.

    Raises
    ------
    ValueError
        If either ``x`` or ``xp`` is not 0 or 1.
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


def single_matrix_element_to_qubit_operator(
    l: int, lp: int, coeff: complex, d: int, encoding: str
) -> QubitOperator:
    """
    Builds a QubitOperator object for a given transition - in other words converts a 
    single bosonic matrix element (coeff * |l><lp|) into a sum of Pauli strings acting on 
    qubits (we need to do this for all transitions to get a big string of Pauli strings!).

    In other words: The total bosonic operator (written in the truncated Fock basis) can be 
    decomposed as a sum over matrix elements |l><lp|. This function takes one such element, 
    encodes each state, and maps the element into qubit space. The result depends on the chosen 
    encoding.

    'single' as h.c. logic in next function.

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

    Notes
    -------
    This is probably the most important function here.

    """
    bits_l = bits_for_level(l, d, encoding)
    bits_lp = bits_for_level(lp, d, encoding)
    support = sorted(bitmask_subset(l, d, encoding) | bitmask_subset(lp, d, encoding)) # all for compact encodings
    local_choices = [] 
    for q in support:
        local_choices.append(one_qubit_map(bits_l[q], bits_lp[q])) # a list of pauli terms (like 1/2( X + iY) for each relevant qubit (in support) 

    op = QubitOperator() 

    for pauli_assignment in product(*local_choices): # * breaks up the list and product expands pauli terms to get all possible strings (each term acts across multiple qubits)
        # pauli_assigment is one possible pauli string (no coefficient yet)
        term = []
        total_coeff = coeff # from position_operator_matrix
        for qubit, (pauli, c) in zip(support, pauli_assignment):
            # zip(support, pauli_assignment): zip({0, 1}, ("X", 0.5)  ("Y", 0.5j)) --> (0, ("X", 0.5)), 1, ("Y", 0.5j)
            # pauli: X/ Y/ Z/ I & c: 0.5 or 0.5j
            total_coeff *= c
            if pauli != "I":
                term.append((qubit, pauli))
        op += QubitOperator(tuple(term), total_coeff)
    return op


def matrix_element_to_qubit_operator(
    l: int, lp: int, coeff: complex, d: int, encoding: str
) -> QubitOperator:
    """
    Map a bosonic matrix element to a Hermitian qubit operator.

    Converts the matrix element

        coeff * |l><lp|

    into its encoded Pauli representation. For off-diagonal elements, the
    Hermitian-conjugate term is added automatically:

        coeff * |l><lp| + coeff.conjugate() * |lp><l|.

    Diagonal elements are included only once. Pauli terms with coefficients
    smaller than ``1e-12`` in magnitude are removed using
    :meth:`QubitOperator.compress`.

    Parameters
    ----------
    l : int
        Bosonic basis-state index associated with the ket ``|l>``.
    lp : int
        Bosonic basis-state index associated with the bra ``<lp|``.
    coeff : complex
        Coefficient multiplying the matrix element ``|l><lp|``.
    d : int
        Dimension of the truncated bosonic Hilbert space.
    encoding : str
        Qubit encoding used to represent the bosonic basis states.

    Returns
    -------
    QubitOperator
        Hermitian Pauli decomposition of the encoded matrix element.

    Notes
    -----
    The phase-space quadrature operators have zero diagonal matrix elements,
    so the diagonal branch is mainly relevant for other bosonic operators
    that may be added to the codebase later.
    """
    op_forward = single_matrix_element_to_qubit_operator(l, lp, coeff, d, encoding)

    if l == lp: # Do not include the hermitian conjugate for the diagonal terms
        op_forward.compress(abs_tol=1e-12)
        return op_forward
    
    op_backward = single_matrix_element_to_qubit_operator(lp, l, coeff.conjugate(), d, encoding)

    op = op_forward + op_backward
    op.compress(abs_tol=1e-12)
    return op


def matrix_to_qubit_operator(mat, d: int, encoding: str, tol: float = 1e-14) -> QubitOperator:
    """
    Convert a full bosonic operator (given as a d x d matrix) into a QubitOperator.

    Parameters
    ----------
    mat : array-like (d x d)
        Matrix representation of the bosonic operator in the truncated Fock basis.
        The bosonic displacement matrix is tri-diagonal (only adjacent states
        connected and so most entries in this are equal to 0).
    d : int
        Bosonic cutoff dimension (d = N_max + 1)
    encoding : str
        Encoding used to map bosonic levels to qubits.
        Supported values: "sb", "gray", "unary".
    tol : float, optional
        Threshold below which matrix elements are treated as zero (default: 1e-14).

    Returns
    -------
    QubitOperator
        Sum of Pauli strings representing the full operator in qubit space.

    Notes
    -------
    
    Only the diagonal and upper-triangular matrix elements are processed.
    For each off-diagonal element, the corresponding lower-triangular
    Hermitian-conjugate term is added automatically by
    ``matrix_element_to_qubit_operator``. Processing the full matrix would
    therefore double count the off-diagonal contributions.
    """
    op = QubitOperator()
    for l in range(d):
        for lp in range(l, d): # only half the triangle
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

    Recall QubitOperator(term, coefficient). 
    The 'term' is a tuple of (qubit_index, Pauli_letter) like ((0, "X"), (2, "Z")) and the 
    'coefficient' is the coefficient in front of this Pauli string. If the 'term' contains
    no entries, it consists entirely on Identities and so it pushed to the end of the 
    "pseudo-alphabetized" Pauli Strings.

    Returns another QubitOperator (just ordered)

    Notes
    -------

    The pseudo-alphabetical ordering used by Sawaya groups structurally similar Pauli strings 
    next to one another, for example terms beginning with X_0, then terms beginning with 
    Y_0, then Z_0, before moving to strings whose first non-identity Pauli appears on a
    higher-index qubit. This ordering is useful because similar neighbouring Pauli strings increase 
    the chance of gate cancellations after the corresponding CNOT ladders are constructed. 
    """

    pauli_priority = {"X": 0, "Y": 1, "Z": 2} # X before Y before Z (for each qubit)

    def key_fn(item):
        term, coeff = item

        if len(term) == 0:
            return (10_000,) # pushes identity terms to the end

        # ensure each Pauli String (multi-qubit operator) is sorted by qubit index
        term = sorted(term, key=lambda x: x[0])

        first_qubit, first_pauli = term[0]
        primary = pauli_priority[first_pauli] # Converts the first Pauli letter into its sorting priority using pauli_priority

        qubit_pattern = tuple(q for q, _ in term)
        # qubit_pattern: [(1, "Z"), (3, "X")] becomes (1, 3)
        pauli_pattern = tuple(pauli_priority[p] for _, p in term) 
        # pauli_pattern: [(1, "Z"), (3, "X")] becomes (2, 0)

        return (primary, qubit_pattern, pauli_pattern)

    # sort terms
    sorted_items = sorted(op.terms.items(), key=key_fn)

    # rebuild operator in that order
    new_op = QubitOperator()
    for term, coeff in sorted_items:
        new_op += QubitOperator(term, coeff)

    return new_op


# op = (
#     QubitOperator(((0, "X"), (1, "Z")), 1.0)
#     + QubitOperator(((0, "Y"), (1, "Z")), 1.0)
#     + QubitOperator(((0, "Z"), (1, "Z")), 1.0)
#     + QubitOperator(((0, "X"), (2, "Z")), 1.0)
#     + QubitOperator(((0, "Y"), (2, "Z")), 1.0)
#     + QubitOperator(((0, "Z"), (2, "Z")), 1.0)
# )

# print("Original operator:")
# for term, coeff in op.terms.items():
#     print(term, coeff)

# print("\nAfter pseudo-alphabetical ordering:")
# sorted_op = pseudo_alphabetical_qubit_operator(op)
# for term, coeff in sorted_op.terms.items():
#     print(term, coeff)
