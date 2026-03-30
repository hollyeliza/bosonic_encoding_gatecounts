import math
import cirq
from openfermion import QubitOperator
from src.pauli_string_formation.mapping import matrix_element_to_qubit_operator, pseudo_alphabetical_qubit_operator


def cirq_qubits(op: QubitOperator) -> dict[int, cirq.Qid]:
    """
    Create a mapping from integer qubit indices appearing in a QubitOperator
    to Cirq LineQubits.

    Parameters
    ----------
    op
        OpenFermion QubitOperator. Its terms are Pauli strings, and each
        Pauli string is a tuple of (qubit_index, pauli_label) pairs.

    Returns
    -------
    dict[int, cirq.Qid]
        Mapping from integer qubit index to the corresponding Cirq qubit.

    Example
    -------
    If op contains:
        X0 Y2 + Z1

    then this returns:
        {
            0: cirq.LineQubit(0),
            1: cirq.LineQubit(1),
            2: cirq.LineQubit(2),
        }
    """

    terms = op.terms # The dictionary inside the QubitOperator that stores each Pauli string and coefficients

    all_indices = set()

    for pauli_string in terms.keys(): 
        # Example pauli_string: ((0, 'Z'), (1, 'X'))
        for local_factor in pauli_string:
            qubit_index = local_factor[0]
            all_indices.add(qubit_index)

    qubit_map = {}
    for i in sorted(all_indices):
        qubit_map[i] = cirq.LineQubit(i)

    return qubit_map

# trial = matrix_element_to_qubit_operator(5, 6, 2, 8, "gray")
# trial_ordered = pseudo_alphabetical_qubit_operator(trial)
# print(f'Pauli string for single transition: {trial_ordered}') # hmm already sorted without pseudo
# trial_2 = cirq_qubits(trial_ordered)
# print(f'cirq mapping: {trial_2}')



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


def cirq_circuit(ordered_terms: QubitOperator) -> tuple[cirq.Circuit, dict[int, cirq.Qid]]:
    """
    Builds a Cirq circuit from an already ordered list of Pauli terms .
    (takes output from sorted_terms_pseudo_alphabetical). It looks through
    all of the Pauli terms in the operator, compiles each term into gates
    and then appends them all to the same circuit.

    It does this by changing basis and then using a CNOT ladder to get the
    parity onto the last qubit.

    Each term is compiled as exp(-i * theta * P), where theta is taken from the 
    coefficient. There should be no non-Hermitian terms because I took the h.c. 
    in mapping.matrix_element_to_qubit_operator.

    Parameters
    ----------
    ordered_terms
        The Pauli string for the bosonic displacement operator in the pseudo - alphabetical
        order

    Returns
    -------
    circuit, qubit_map
    """

    qubit_map = cirq_qubits(ordered_terms) # Creates all the required cirq qubits
    circuit = cirq.Circuit()

    for term, coeff in ordered_terms.terms.items():
        # term = ((0, 'X'), (2, 'Y'))
        # coeff = 0.3
        theta = coeff.real # the phase angle used for the globl phase on last qubit after all entangled with CNOTs

        qubits = [qubit_map[q] for q, _ in term] # the cirq qubit associted with the term
        paulis = [p for _, p in term] # the paulis associated with the cirq qubits in this term

        # 1. Basis change → Z basis
        for qubit, pauli in zip(qubits, paulis):
            basis_change_into_z(circuit, qubit, pauli)

        # 2–3. Parity + Z rotation (CNOT ladder)
        if len(qubits) == 1:
            circuit.append(cirq.rz(2 * theta).on(qubits[0]))
        else:
            # Forward ladder: q0->q1, q1->q2, ..., q_{n-2}->q_{n-1}
            for control, target in zip(qubits[:-1], qubits[1:]):
                circuit.append(cirq.CNOT(control, target))

            # Apply the phase on the last qubit
            circuit.append(cirq.rz(2 * theta).on(qubits[-1]))

            # Undo the ladder in reverse
            for control, target in reversed(list(zip(qubits[:-1], qubits[1:]))):
                circuit.append(cirq.CNOT(control, target))

        # 4. Undo basis change
        for qubit, pauli in reversed(list(zip(qubits, paulis))):
            undo_basis_change(circuit, qubit, pauli)

    return circuit, qubit_map

# trial_3_circuit,  trial_3_qubit_map = cirq_circuit(trial) # Trial for one transition
# print(type(trial_3_circuit))


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




# circuit is being treated as QubitOoperator but circuit is a cirq.circuit with these attributes...

# print(f'The number of cnots in the trial before cirq optimization: {count_cnots(trial_3_circuit)}')

# optimized_circuit_trial = optimize_cirq_circuit(trial_3_circuit)
# print(f'The number of cnots in the trial after cirq optimization: {count_cnots(optimized_circuit_trial)}')
