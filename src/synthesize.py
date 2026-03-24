import numpy as np
from qiskit import QuantumCircuit, transpile


def apply_basis_change(qc: QuantumCircuit, qubit: int, pauli: str):
    if pauli == "X":
        qc.h(qubit)
    elif pauli == "Y":
        qc.sdg(qubit)
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
    Append exp(-i theta P) for one Pauli term P using a CNOT ladder.
    Global phase/sign conventions are not important for gate counting here.
    """
    if len(term) == 0:
        return

    qubits = [q for q, p in term]
    paulis = [p for q, p in term]

    for q, p in zip(qubits, paulis):
        apply_basis_change(qc, q, p)

    if len(qubits) == 1:
        qc.rz(2.0 * theta, qubits[0])
    else:
        target = qubits[-1]
        for control in qubits[:-1]:
            qc.cx(control, target)
        qc.rz(2.0 * theta, target)
        for control in reversed(qubits[:-1]):
            qc.cx(control, target)

    for q, p in reversed(list(zip(qubits, paulis))):
        undo_basis_change(qc, q, p)


def cnot_count_of_compiled_circuit(qc: QuantumCircuit, optimization_level: int = 3) -> int:
    tqc = transpile(qc, basis_gates=["cx", "rz", "h", "s", "sdg"], optimization_level=optimization_level)
    counts = tqc.count_ops()
    return int(counts.get("cx", 0)), tqc