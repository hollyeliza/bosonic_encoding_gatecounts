# Bosonic Encoding Gate Counts

Reproduction of CNOT gate count scaling for the bosonic position operator
across different encodings (standard binary, Gray, unary), following:

Sawaya et al., npj Quantum Information (2020)

The unoptimized CNOT counts are looked at alongside different optimization methods
(the way the paper described and also qiskit and cirq) to try and recreate the CNOT
counts obtained in the paper for the encodings.

## What this does

- Constructs truncated bosonic displacement operator
- Maps it to qubit operators using different encodings
- Expands into Pauli strings using OpenFermion
- Synthesizes circuits and counts CNOT gates (naive equation way and using compilation)
- 'Tries to' reproduce scaling trends in the paper

## File structure

Note: Since all changed - need to update

src/
- encodings_b.py        # bit encodings (unary, Gray, binary)
- mapping.py            # maps |l><l'| → Pauli strings
- bosonic_ops.py        # constructs bosonic operators (e.g. q)
- gate_count.py         # CNOT counting utilities

tests/
- test_mapping.py       # unit tests for mapping functions


## Key idea

Bosonic operators are sparse in the Fock basis but become non-local 
when mapped to qubits. This project quantifies the resulting gate overhead.

Mappings (mapping.py)

We encode the bosonic state (phonon occupation) from |l> to |R(l)> where
R(l) represents the encoding (gray, sb, unary). To apply the bosonic
operators like the displacement operator, we need |R(l)><R(l')|.
To get this we break this outer product into a tensor product of one
qubit outer products, that each have a nice mpping to a Pauli gate
As the qubit can be 0 or 1, there are only 4 possibilities

Example:

|5><6|
#Gray: |0111><0101| Note: sequential occupations differ by 1 qubit - local!
0111><0101| = |0><0| x|1><1| x |1><0| x |1><1|

Each of these one qubit outer products can be replaced in in terms of X,Y,Z
#operators and this forms a sum of Pauli strings

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

## Example

Construct the displacement operator and map to qubits:

q = bosonic_disp_operator_matrix(d=8)
op = matrix_to_qubit_operator(q, d=8, encoding="gray")
print(op)


## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
