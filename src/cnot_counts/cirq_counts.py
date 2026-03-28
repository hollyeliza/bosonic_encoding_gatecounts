import math
import cirq
from openfermion import QubitOperator
from optimize.cirq_comp_op import qubit_operator_to_trotter_circuit, optimize_cirq_circuit, count_cnots

def cirq_cnot_count_after_optimization(
    op: QubitOperator,
    time: float = 1.0,
    term_order: str = "default",
) -> tuple[int, int, cirq.Circuit, cirq.Circuit]:
    """
    Build the Cirq Trotter circuit, optimize it, and return:

        raw_cnot_count, optimized_cnot_count, raw_circuit, optimized_circuit
    """
    raw_circuit, qubit_map = qubit_operator_to_trotter_circuit(
        op,
        time=time,
        term_order=term_order,
    )

    optimized_circuit = optimize_cirq_circuit(raw_circuit)

    raw_count = count_cnots(raw_circuit)
    optimized_count = count_cnots(optimized_circuit)

    return raw_count, optimized_count, raw_circuit, optimized_circuit