# The optimize files contain functions that take the Pauli strings for the entire 
# displacement operator (in the correct order - pseudoalphabetical), compile to 
# a circuit and optimizes this circuit, cutting down the number of required
# CNOTs

import numpy as np
from qiskit import QuantumCircuit, transpile


def apply_basis_change(qc: QuantumCircuit, qubit: int, pauli: str):
    if pauli == "X":
        qc.h(qubit)
    elif pauli == "Y":
        qc.sdg(qubit) # s gate with phase - undoes s
        qc.h(qubit)
    elif pauli == "Z":
        pass
    else:
        raise ValueError(f"Unexpected Pauli: {pauli}")


def undo_basis_change(qc: QuantumCircuit, qubit: int, pauli: str):
    if pauli == "X":
        qc.h(qubit)
    elif pauli == "Y":
        qc.h(qubit)
        qc.s(qubit)
    elif pauli == "Z":
        pass
    else:
        raise ValueError(f"Unexpected Pauli: {pauli}")


def qiskit_circuit(
    qc: QuantumCircuit,
    term,
    theta: float,
) -> None:
    """
    Append the unitary exp(-i * theta * P) for a single Pauli string P
    to an existing Qiskit circuit.

    Here P is a product of Pauli operators acting on different qubits,
    for example X_0 Y_1 Z_2. This corresponds to one term in a Hamiltonian
    after Pauli decomposition.

    The implementation proceeds by:
    1. Rotating each involved qubit into the Z basis,
    2. Using a CNOT ladder to accumulate the parity of the Pauli string
       onto the last qubit,
    3. Applying a single Rz(2 * theta) rotation on that qubit,
    4. Uncomputing the CNOT ladder,
    5. Undoing the basis changes.

    Parameters
    ----------
    qc
        QuantumCircuit to append to.
    term
        Pauli-string term, e.g. ((0, 'X'), (1, 'Y'), (2, 'Z')).
    theta
        Rotation angle in exp(-i * theta * P).

    Returns
    -------
    None
        The circuit is modified in place.

    Notes
    -----
    Identity terms produce only a global phase and are skipped here.
    Global phase/sign conventions are ignored, as they do not affect
    gate counting.
    """
    if len(term) == 0:
        return

    qubits = [q for q, _ in term]
    paulis = [p for _, p in term]

    # 1. Rotate into the Z basis.
    for q, p in zip(qubits, paulis):
        apply_basis_change(qc, q, p)

    # 2–4. CNOT ladder + Z rotation + uncompute.
    if len(qubits) == 1:
        qc.rz(2.0 * theta, qubits[0])
    else:
        # Forward ladder: q0->q1, q1->q2, ..., q_{n-2}->q_{n-1}
        for control, target in zip(qubits[:-1], qubits[1:]):
            qc.cx(control, target)

        # Apply phase on the last qubit
        qc.rz(2.0 * theta, qubits[-1])

        # Undo ladder in reverse
        for control, target in reversed(list(zip(qubits[:-1], qubits[1:]))):
            qc.cx(control, target)

    # 5. Undo basis changes.
    for q, p in reversed(list(zip(qubits, paulis))):
        undo_basis_change(qc, q, p)


# Question: Is the hermitian conjugate counted in the cost - this could really affect things!!
# Also explicitly write a function to calculate the naive count also.

def transpilation_counts(qc: QuantumCircuit, optimization_level: int = 3) -> int:
    """
    Transpiles the quantum circuit and returns the number of CNOT (CX) gates after 
    optimization.

    The circuit is compiled using Qiskit's transpiler into a basis consisting
    of {"cx", "rz", "h", "s", "sdg"}, which reflects a typical hardware-native
    gate set for superconducting qubits.

    The optimization level controls how aggressively Qiskit simplifies the circuit:
        - 0: no optimization
        - 1–3: increasing levels of gate cancellation, merging, and re-synthesis

    This function is used to obtain a more realistic estimate of the entangling
    gate cost (CNOT count) compared to naive analytical estimates such as
    2(p-1) per Pauli string.
    """
    tqc = transpile(qc, basis_gates=["cx", "rz", "h", "s", "sdg"], optimization_level=optimization_level)
    counts = tqc.count_ops()
    return int(counts.get("cx", 0)), tqc