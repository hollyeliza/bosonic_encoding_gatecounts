# This is currently specific to qiskit
# Make specific - make sure you sort alphabetical and prep according to SI before passing
# into the optimizer

from qiskit import QuantumCircuit
from openfermion import QubitOperator

from src.optimize.qiskit_comp_op import qiskit_circuit, optimize_qiskit_circuit


def count_cnots_qiskit(circuit: QuantumCircuit) -> int:
    """
    Count CNOT gates in a Qiskit circuit.
    """
    total = 0
    for instruction, qargs, cargs in circuit.data:
        if instruction.name == "cx":
            total += 1
    return total


def qiskit_cnot_count_before_and_after_optimization(
    op: QubitOperator,
) -> tuple[int, int]:
    """
    Build the Qiskit circuit, optimize it, and return:

        raw_cnot_count, optimized_cnot_count
    """
    pre_qiskit_op_circuit = qiskit_circuit(op)
    pre_op_count = count_cnots_qiskit(pre_qiskit_op_circuit)

    optimized_count, qiskit_op_circuit = optimize_qiskit_circuit(pre_qiskit_op_circuit)

    return pre_op_count, optimized_count