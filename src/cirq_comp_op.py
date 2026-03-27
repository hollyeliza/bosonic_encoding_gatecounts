# I haven't yet compared qiskit and cirq compilation methods

import math
import cirq
from openfermion import QubitOperator


def qubits_from_qubit_operator(op: QubitOperator) -> dict[int, cirq.Qid]:
    """
    Create a consistent mapping from OpenFermion qubit indices to Cirq qubits.
    """
    all_indices = sorted({q for term in op.terms for q, _ in term})
    return {i: cirq.LineQubit(i) for i in all_indices}


def basis_change_into_z(circuit: cirq.Circuit, qubit: cirq.Qid, pauli: str) -> None:
    """
    Append basis-change gates so that measuring/applying a Z-phase in the new basis
    is equivalent to acting in the original Pauli basis.

    X -> H Z H
    Y -> S^-1 H Z H S
    Z -> already Z basis
    """
    if pauli == "X":
        circuit.append(cirq.H(qubit))
    elif pauli == "Y":
        circuit.append(cirq.S(qubit) ** -1)
        circuit.append(cirq.H(qubit))
    elif pauli == "Z":
        pass
    else:
        raise ValueError(f"Unsupported Pauli: {pauli}")


def undo_basis_change(circuit: cirq.Circuit, qubit: cirq.Qid, pauli: str) -> None:
    """
    Undo the basis change applied in basis_change_into_z.
    """
    if pauli == "X":
        circuit.append(cirq.H(qubit))
    elif pauli == "Y":
        circuit.append(cirq.H(qubit))
        circuit.append(cirq.S(qubit))
    elif pauli == "Z":
        pass
    else:
        raise ValueError(f"Unsupported Pauli: {pauli}")


def pauli_term_to_cirq_circuit(
    term: tuple[tuple[int, str], ...],
    angle: float,
    qubit_map: dict[int, cirq.Qid],
) -> cirq.Circuit:
    """
    Compile exp(-i * angle * P) for a single Pauli string P into a Cirq circuit.

    Parameters
    ----------
    term
        OpenFermion term representation, e.g. ((0, 'X'), (2, 'Y')).
    angle
        The angle in exp(-i * angle * P).
        Since Cirq's rz(phi) implements exp(-i phi Z / 2), we apply rz(2*angle).
    qubit_map
        Mapping from integer indices to Cirq qubits.

    Returns
    -------
    cirq.Circuit
        Circuit implementing the Pauli-string rotation.
    """
    circuit = cirq.Circuit()

    if len(term) == 0:
        # Identity term: global phase only, no CNOT cost.
        return circuit

    qubits = [qubit_map[q] for q, _ in term]
    paulis = [p for _, p in term]

    # 1. Change basis so every local Pauli becomes Z.
    for qubit, pauli in zip(qubits, paulis):
        basis_change_into_z(circuit, qubit, pauli)

    # 2. CNOT ladder to collect parity onto last qubit.
    for control, target in zip(qubits[:-1], qubits[1:]):
        circuit.append(cirq.CNOT(control, target))

    # 3. Z rotation on the last qubit.
    # exp(-i * angle * Z) = Rz(2*angle)
    circuit.append(cirq.rz(2 * angle).on(qubits[-1]))

    # 4. Undo CNOT ladder.
    for control, target in reversed(list(zip(qubits[:-1], qubits[1:]))):
        circuit.append(cirq.CNOT(control, target))

    # 5. Undo basis changes.
    for qubit, pauli in reversed(list(zip(qubits, paulis))):
        undo_basis_change(circuit, qubit, pauli)

    return circuit


def qubit_operator_to_trotter_circuit(
    op: QubitOperator,
    time: float = 1.0,
    term_order: str = "default",
) -> tuple[cirq.Circuit, dict[int, cirq.Qid]]:
    """
    Build a first-order Trotter circuit for exp(-i * time * H), where
    H is a Hermitian OpenFermion QubitOperator with real coefficients.

    Each term c*P is compiled as exp(-i * time * c * P).

    Parameters
    ----------
    op
        Hermitian QubitOperator with real coefficients.
    time
        Evolution time.
    term_order
        "default" keeps OpenFermion iteration order.
        "pseudo_alpha" sorts by qubit index then Pauli type X<Y<Z,
        roughly mimicking the ordering used in some gate-count papers.

    Returns
    -------
    circuit, qubit_map
    """
    qubit_map = qubits_from_qubit_operator(op)

    items = list(op.terms.items())

    if term_order == "pseudo_alpha":
        pauli_priority = {"X": 0, "Y": 1, "Z": 2}

        def sort_key(item):
            term, coeff = item
            if len(term) == 0:
                return ((10_000, 10_000),)
            return tuple((q, pauli_priority[p]) for q, p in term)

        items = sorted(items, key=sort_key)

    circuit = cirq.Circuit()

    for term, coeff in items:
        # Identity term only contributes a global phase; ignore for CNOT counts.
        if len(term) == 0:
            continue

        if abs(coeff.imag) > 1e-12:
            raise ValueError(
                f"Term {term} has complex coefficient {coeff}. "
                "For circuit compilation as time evolution, use a Hermitian operator "
                "with real coefficients."
            )

        angle = time * coeff.real
        circuit += pauli_term_to_cirq_circuit(term, angle, qubit_map)

    return circuit, qubit_map


def count_cnots(circuit: cirq.Circuit) -> int:
    """
    Count CNOT gates in a Cirq circuit.
    """
    return sum(
        1
        for op in circuit.all_operations()
        if isinstance(op.gate, cirq.CNotPowGate) and abs(op.gate.exponent - 1) < 1e-12
    )


def optimize_cirq_circuit(circuit: cirq.Circuit) -> cirq.Circuit:
    """
    Apply a few simple Cirq optimizations. This is not guaranteed to reproduce
    a paper's custom optimizer exactly, but it is a reasonable baseline.
    """
    optimized = circuit.copy()

    # Merge / simplify 1q gates and eject Z's where possible.
    cirq.merge_k_qubit_unitaries(optimized, k=1, rewriter=lambda op: op)
    cirq.eject_z(optimized)
    cirq.drop_negligible_operations(optimized)
    cirq.drop_empty_moments(optimized)

    return optimized