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

    These operations rotate into the Z basis, we do the operations, and then undo.
    This function takes us from X to Z (requires 1 H) or Y to Z (requires 1 S^-1 and 1 H)

    (from here) X -> H (to here) Z H
    (from here) Y -> S^-1 H (to here) Z H S
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

    These operations now take us out of the Z basis.
    This function takes us from Z back to X (requires 1 H) or Z to Y (requires 1 H and 1 S).

    X -> H (from here) Z H (to here)
    Y -> S^-1 H (from here) Z H S (to here)
    Z -> already Z basis 
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


def cirq_circuit(
    ordered_terms,
) -> tuple[cirq.Circuit, dict[int, cirq.Qid]]:
    """
    Build a Cirq circuit from an already ordered list of Pauli terms.

    Each term is compiled as exp(-i * theta * P), where theta is taken
    directly from the coefficient.

    Parameters
    ----------
    ordered_terms
        Iterable of (term, coeff) pairs, already in the desired order.
        The coefficient is interpreted directly as the rotation angle.

    Returns
    -------
    circuit, qubit_map
    """
    all_indices = sorted({q for term, _ in ordered_terms for q, _ in term})
    qubit_map = {i: cirq.LineQubit(i) for i in all_indices}

    circuit = cirq.Circuit()

    for term, coeff in ordered_terms:
        if len(term) == 0:
            continue

        if abs(coeff.imag) > 1e-12:
            raise ValueError(
                f"Term {term} has complex coefficient {coeff}. "
                "Expected real coefficients for unitary evolution."
            )

        theta = coeff.real

        qubits = [qubit_map[q] for q, _ in term]
        paulis = [p for _, p in term]

        # 1. Basis change → Z basis
        for qubit, pauli in zip(qubits, paulis):
            basis_change_into_z(circuit, qubit, pauli)

        # 2–3. Parity + Z rotation (star pattern)
        if len(qubits) == 1:
            circuit.append(cirq.rz(2 * theta).on(qubits[0]))
        else:
            target = qubits[-1]

            for control in qubits[:-1]:
                circuit.append(cirq.CNOT(control, target))

            circuit.append(cirq.rz(2 * theta).on(target))

            for control in reversed(qubits[:-1]):
                circuit.append(cirq.CNOT(control, target))

        # 4. Undo basis change
        for qubit, pauli in reversed(list(zip(qubits, paulis))):
            undo_basis_change(circuit, qubit, pauli)

    return circuit, qubit_map


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

def count_cnots(circuit: cirq.Circuit) -> int:
    """
    Count exact CNOT gates in a Cirq circuit.
    """
    total = 0
    for op in circuit.all_operations():
        if isinstance(op.gate, cirq.CNotPowGate) and abs(op.gate.exponent - 1) < 1e-12:
            total += 1
    return total


def optimize_cirq_circuit(circuit: cirq.Circuit) -> cirq.Circuit:
    """
    Apply a reasonable set of Cirq optimizations and return the optimized circuit.

    This is a Cirq baseline optimizer, not the same as the paper's custom optimizer.
    """
    optimized = circuit.copy()

    # Run several passes because some simplifications expose others.
    for _ in range(3):
        cirq.expand_composite(optimized)

        # Clean up single-qubit structure.
        cirq.merge_single_qubit_gates_to_phxz(optimized)
        cirq.eject_z(optimized)

        # Remove tiny / empty leftovers.
        cirq.drop_negligible_operations(optimized)
        cirq.drop_empty_moments(optimized)

        # Try to merge adjacent unitary blocks on 1 and 2 qubits.
        cirq.merge_k_qubit_unitaries(
            optimized,
            k=1,
            rewriter=lambda op: op
        )
        cirq.merge_k_qubit_unitaries(
            optimized,
            k=2,
            rewriter=lambda op: op
        )

        cirq.drop_negligible_operations(optimized)
        cirq.drop_empty_moments(optimized)

    return optimized


