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


def append_pauli_evolution(
    qc: QuantumCircuit,
    term,
    theta: float,
):
    """
    Append the unitary exp(-i θ P) for a single Pauli string P (e.g. X₀Y₁Z₂...).

    Here P is a *product* of Pauli operators acting on different qubits, not a sum.
    This corresponds to one term in the Hamiltonian after decomposition.

    The theta is from the combo of coefficients (squre roots) and 1/2's in front of Pauli 
    operator for single qubit transition).

    The implementation proceeds by:
    1. Rotating each qubit into the Z basis (so P → Z⊗Z⊗...),
    2. Using a CNOT ladder to compute the parity of the involved qubits onto a target,
    3. Applying a single Rz rotation on the target qubit, which encodes the phase
    exp(-i θ P) for the entire string,
    4. Uncomputing the parity and undoing the basis change.

    Global phase/sign conventions are ignored, as they do not affect gate counting.
    """
    if len(term) == 0:
        return

    qubits = [q for q, p in term] # Just gets the qubit number the paulis in the term correspond to
    paulis = [p for q, p in term]

    for q, p in zip(qubits, paulis):
        apply_basis_change(qc, q, p) # We can only implement Z rotation

    if len(qubits) == 1:
        qc.rz(2.0 * theta, qubits[0]) # No CNOT cascade required! Just apply the rotation.
    else: # Okay so CNOT cascade required
        target = qubits[-1] # Last qubit accumulates all of the phase
        for control in qubits[:-1]:
            qc.cx(control, target) # Adds all the qubits - assumes full connectivity - trapped ion
        qc.rz(2.0 * theta, target) # Z rotation in single qubit
        for control in reversed(qubits[:-1]):
            qc.cx(control, target) #  Uncomputes the parity

    for q, p in reversed(list(zip(qubits, paulis))):
        undo_basis_change(qc, q, p) # undo the basis change


def cnot_count_of_compiled_circuit(qc: QuantumCircuit, optimization_level: int = 3) -> int:
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