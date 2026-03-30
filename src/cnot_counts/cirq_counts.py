import math
import cirq
from openfermion import QubitOperator
from src.optimize.cirq_comp_op import cirq_circuit, optimize_cirq_circuit


def count_cnots_cirq(circuit: cirq.Circuit) -> int:
    """
    Count exact CNOT gates in a Cirq circuit.
    Note: cirq circuit returns a tuple where the first term is the actual cirq circuit. Make sure
    that you feed the first term only into this function.
    """
    total = 0
    for op in circuit.all_operations():
        if isinstance(op.gate, cirq.CNotPowGate) and abs(op.gate.exponent - 1) < 1e-12:
            total += 1
    return total


def cirq_cnot_count_before_and_after_optimization(op: QubitOperator) -> tuple[int, int, cirq.Circuit, cirq.Circuit]:
    """
    Build the Cirq Trotter circuit, optimize it, and return:

        raw_cnot_count, optimized_cnot_count, raw_circuit, optimized_circuit
    """

    # Remember that cirq_circuit returns a tuple so we need to first unpack the 
    # cirq circuit and qubit_map. Although qubit map is not really used so maybe
    # it is better for cirq_circuit to just return the actual circuit lol

    pre_cirq_op_circuit, cirq_qubit_map = cirq_circuit(op)

    cirq_op_circuit = optimize_cirq_circuit(pre_cirq_op_circuit)

    pre_op_count = count_cnots_cirq(pre_cirq_op_circuit)
    optimized_count = count_cnots_cirq(cirq_op_circuit)

    return pre_op_count, optimized_count