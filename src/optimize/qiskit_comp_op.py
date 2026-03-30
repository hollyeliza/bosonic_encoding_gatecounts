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


from qiskit import QuantumCircuit
from openfermion import QubitOperator


def append_qiskit_pauli_term(
    qc: QuantumCircuit,
    term,
    theta: float,
) -> None:
    """
    Append exp(-i * theta * P) for a single Pauli string P
    to an existing Qiskit circuit.

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
    """
    if len(term) == 0:
        return

    qubits = [q for q, _ in term]
    paulis = [p for _, p in term]

    # 1. Rotate into the Z basis
    for q, p in zip(qubits, paulis):
        apply_basis_change(qc, q, p)

    # 2–4. CNOT ladder + Z rotation + uncompute
    if len(qubits) == 1:
        qc.rz(2.0 * theta, qubits[0])
    else:
        for control, target in zip(qubits[:-1], qubits[1:]):
            qc.cx(control, target)

        qc.rz(2.0 * theta, qubits[-1])

        for control, target in reversed(list(zip(qubits[:-1], qubits[1:]))):
            qc.cx(control, target)

    # 5. Undo basis changes
    for q, p in reversed(list(zip(qubits, paulis))):
        undo_basis_change(qc, q, p)


def qiskit_circuit(ordered_terms: QubitOperator) -> QuantumCircuit:
    """
    Build a full Qiskit circuit from a QubitOperator whose terms are already
    in the desired iteration order.

    Each term c * P is compiled as exp(-i * theta * P), where theta is taken
    from the real part of the coefficient c.

    Parameters
    ----------
    ordered_terms
        QubitOperator whose terms are already ordered as desired
        (for example via pseudo_alphabetical_qubit_operator).

    Returns
    -------
    QuantumCircuit
        The full Trotter-style circuit obtained by appending all term circuits.
    """
    # Find how many qubits are needed
    max_qubit = -1
    for term in ordered_terms.terms:
        for q, _ in term:
            max_qubit = max(max_qubit, q)

    n_qubits = max_qubit + 1 if max_qubit >= 0 else 1
    qc = QuantumCircuit(n_qubits)

    for term, coeff in ordered_terms.terms.items():
        theta = coeff.real
        append_qiskit_pauli_term(qc, term, theta)

    return qc


# Question: Is the hermitian conjugate counted in the cost - this could really affect things!!
# Also explicitly write a function to calculate the naive count also.

from qiskit import QuantumCircuit, transpile


def optimize_qiskit_circuit(
    qc: QuantumCircuit,
    optimization_level: int = 3,
) -> tuple[int, QuantumCircuit]:
    """
    Transpile the quantum circuit and return:

        optimized_cnot_count, optimized_circuit

    The circuit is compiled into the basis
    {"cx", "rz", "h", "s", "sdg"}.

    Parameters
    ----------
    qc
        Input circuit.
    optimization_level
        Qiskit transpiler optimization level:
            0: no optimization
            1-3: increasing optimization effort

    Returns
    -------
    optimized_cnot_count, optimized_circuit
    """
    tqc = transpile(
        qc,
        basis_gates=["cx", "rz", "h", "s", "sdg"],
        optimization_level=optimization_level,
    )
    counts = tqc.count_ops()
    return int(counts.get("cx", 0)), tqc