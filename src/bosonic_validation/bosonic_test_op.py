from __future__ import annotations

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

from src.pauli_string_formation.mapping import (
    matrix_to_qubit_operator,
    pseudo_alphabetical_qubit_operator,
)
from src.pauli_string_formation.ps_quadratures import position_operator_matrix
from src.optimize.qiskit_comp_op import (qiskit_circuit, qiskit_optimizer_and_counts)

import numpy as np

def position_operator_test_circuit(
    d: int = 16,
    encoding: str = "sb",
    theta: float = 1.0,
    time: float = 1.0,
    trotter_steps: int = 1,
    optimize: bool = False,
) -> QuantumCircuit:
    """
    Construct a qubit circuit approximating

        U(t) = exp(-i theta G t),

    where G is the truncated bosonic position operator.

    The operator is mapped to a sum of Pauli strings and implemented using
    a first-order Trotter product formula.

    Parameters
    ----------
    d
        Dimension of the truncated bosonic Fock space. For standard binary
        or Gray encoding, the required number of qubits is approximately
        log2(d). For unary encoding, d qubits are used.

    encoding
        Bosonic-to-qubit encoding, such as ``"sb"``, ``"gray"``, or
        ``"unary"``.

    theta
        Coupling strength multiplying the position operator (leave as 1)

    time
        Total evolution time.

    trotter_steps
        Number of first-order Trotter steps.

    optimize
        Whether to optimize the resulting circuit using the Qiskit
        transpiler (leave as true)

    Returns
    -------
    QuantumCircuit
        Circuit approximating exp(-i theta G t).

    """

    pos_op_matrix = position_operator_matrix(d)
    qubit_operator = matrix_to_qubit_operator(pos_op_matrix, d, encoding)

    # Scale the generator by theta (leaving as 1 in the tests so not that important)
    qubit_operator *= theta

    ordered_operator = pseudo_alphabetical_qubit_operator(qubit_operator)

    circuit = qiskit_circuit(ordered_operator, time=time, trotter_steps=trotter_steps)

    if optimize:
        _, circuit = qiskit_optimizer_and_counts(circuit)

    return circuit


# The following 3 functions are decoder functions to check that the different encodings are being implemeted...

# Why is a decoder needed you may ask

# In bosonic qiskit, the wigner functions expect states in standard binary and so after evolving the state in the 
# chosen encoding, we need to apply a decoding circuit, so that the same state is represented in standard binary and 
# can be plot on the wigner functions


def gray_to_binary_circuit(num_qubits: int) -> QuantumCircuit:
    """
    This circuit is applied after evolving the state with a Gray-encoded
    bosonic circuit. It rearranges the computational-basis states so that
    each amplitude appears at the binary index corresponding to the same
    Fock state, as required for plotting the Wigner function.

    For example, the Gray-code state |011> representing Fock state |2>
    is mapped to the binary state |010>. The amplitude itself is unchanged.

    Qubit 0 is assumed to be the least-significant bit.

    If the higher bit is 1, flip the bit below it. If not leave it. This
    converts to binary/

    Parameters
    ----------
    num_qubits
        Number of qubits in the encoded register.

    Returns
    -------
    QuantumCircuit
        The Gray-to-binary decoder circuit.
    """
    qc = QuantumCircuit(num_qubits) # empty circuit with correct number of bits

    for control in range(num_qubits - 1, 0, -1):
        qc.cx(control, control - 1)

    return qc


# The same can't be done for unary because of the different number of qubits
# Instead you have to take the amplitudes from each Fock state and make a new 
# state vector for plotting


def decode_unary_state(
    custom_circuit: QuantumCircuit,
    d: int,
) -> Statevector:
    """
    Simulate a unary-encoded bosonic circuit and convert the result into
    a standard Fock-basis statevector for plotting the Wigner function.

    In unary encoding, Fock state |n> is represented by the computational
    basis state in which only qubit n is equal to 1. The circuit is applied
    to the unary vacuum, and the amplitudes of the valid unary states are
    extracted in Fock-state order:

        [c_0, c_1, ..., c_{d-1}].

    Parameters
    ----------
    custom_circuit
        Circuit acting on the unary-encoded bosonic register.

    d
        Dimension of the truncated Fock space. Unary encoding uses d qubits.

    Returns
    -------
    Statevector
        Length-d statevector containing the amplitudes of the Fock states
        |0>, |1>, ..., |d-1>.
    """
    if custom_circuit.num_qubits != d:
        raise ValueError(
            f"For unary encoding, cutoff d={d} requires {d} qubits, "
            f"but the circuit uses {custom_circuit.num_qubits}."
        )

    # Basis indices of the valid one-hot unary states.
    unary_indices = [1 << n for n in range(d)]

    # Prepare the unary vacuum:
    # |0>_Fock is represented by qubit 0 being equal to 1.
    initial_state = np.zeros(2**d, dtype=complex)
    initial_state[unary_indices[0]] = 1.0

    # Evolve the unary vacuum.
    unary_state = Statevector(initial_state).evolve(custom_circuit)

    # Extract the amplitudes of |0>, |1>, ..., |d-1>.
    fock_amplitudes = np.array(
        [unary_state.data[index] for index in unary_indices],
        dtype=complex,
    )

    return Statevector(fock_amplitudes)
