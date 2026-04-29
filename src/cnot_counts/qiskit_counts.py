# This is currently specific to qiskit
# Make specific - make sure you sort alphabetical and prep according to SI before passing
# into the optimizer

from qiskit import QuantumCircuit
from openfermion import QubitOperator
from src.optimize.qiskit_comp_op import qiskit_circuit, qiskit_optimizer_and_counts

def qiskit_cnot_count_before_and_after_optimization(
    op: QubitOperator,
) -> tuple[int, int]:
    """
    Build the Qiskit circuit for a QubitOperator and return the CNOT count
    before and after Qiskit transpiler optimization.
    """
    circuit = qiskit_circuit(op)
    pre_count = int(circuit.count_ops().get("cx", 0))

    optimized_count, _ = qiskit_optimizer_and_counts(circuit)

    return pre_count, optimized_count