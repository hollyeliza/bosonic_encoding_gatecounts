# This is currently specific to qiskit
# Make specific - make sure you sort alphabetical and prep according to SI before passing
# into the optimizer

import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit

from pauli_string_formation.bosonic_disp import bosonic_disp_operator_matrix
from pauli_string_formation.mapping import (
    matrix_to_qubit_operator,
    naive_cnot_count_from_qubit_operator,
    sorted_terms_pseudo_alphabetical,
)
from pauli_string_formation.encodings_b import n_qubits
from optimize.qiskit_comp_op import append_pauli_evolution, cnot_count_of_compiled_circuit

def build_full_circuit_from_operator(op, nq: int, theta_scale: float = 1.0):
    """
    Construct a circuit for a Pauli-operator decomposition by appending
    exp(-i θ_k P_k) for each Pauli term P_k in a fixed order.

    Here `op` is a QubitOperator of the form
        op = sum_k coeff_k * P_k.

    For each term:
        θ_k = coeff_k * theta_scale

    The resulting circuit is intended for gate-count estimation rather than
    high-accuracy simulation, so only the real part of each coefficient is used.

    Warning: Not sure about QubitOperator implementation.
    """
    qc = QuantumCircuit(nq)
    for term, coeff in sorted_terms_pseudo_alphabetical(op):
        if abs(coeff) < 1e-12:
            continue
        theta = float(np.real(coeff)) * theta_scale # This sets the evolution angle for that Pauli term.
        # For q-hat all coefficients should be real after Hermitian summation.
        append_pauli_evolution(qc, term, theta)
    return qc


def run_counts(d_values):
    """
    For each cutoff d in d_values and for each encoding ("sb", "gray", "unary"),
    compute:

    1. the encoded Pauli representation of the bosonic position operator q,
    2. the naive CNOT count using the rule 2(p-1) for each Pauli string,
    3. the compiled CNOT count obtained by building and transpiling a circuit.

    Returns
    -------
    naive : dict[str, list[int]]
        Naive CNOT counts for each encoding as a function of d.
    compiled : dict[str, list[int]]
        Transpiled CNOT counts for each encoding as a function of d.
    """
    encodings = ["sb", "gray", "unary"]

    compiled = {enc: [] for enc in encodings} # using qiskit

    for d in d_values:
        print(f"Running d={d}")
        qmat = bosonic_disp_operator_matrix(d)

        for enc in encodings:
            op = matrix_to_qubit_operator(qmat, d=d, encoding=enc)
            nq = n_qubits(d, enc)
            qc = build_full_circuit_from_operator(op, nq=nq, theta_scale=1.0)
            compiled_count = cnot_count_of_compiled_circuit(qc, optimization_level=3)

            compiled[enc].append(compiled_count)

            print(
                f"  {enc:5s} | nq={nq:2d} | #terms={len(op.terms):4d} | "
                f"compiled CX={compiled_count:5d}"
            )
    return compiled