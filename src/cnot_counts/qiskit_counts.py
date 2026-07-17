from openfermion import QubitOperator
from src.optimize.qiskit_comp_op import qiskit_circuit, qiskit_optimizer_and_counts

def qiskit_cnot_count_before_and_after_op(
    op: QubitOperator,
    time: float = 1.0,
    trotter_steps: int = 1,
) -> tuple[int, int]:
    """
    Build the Qiskit circuit for a QubitOperator and return the CNOT count
    before and after Qiskit transpiler optimization.

    Parameters
    ----------
    op
        QubitOperator to compile.

    time
        Total evolution time.

    trotter_steps
        Number of first-order Trotter steps used to approximate the
        time-evolution operator.

    Returns
    -------
    tuple[int, int]
        The CNOT count before optimization and after Qiskit transpiler
        optimization.
    """
    circuit = qiskit_circuit(op, time=time, trotter_steps=trotter_steps)

    pre_op_count = int(circuit.count_ops().get("cx", 0))

    optimized_count, _ = qiskit_optimizer_and_counts(circuit)

    return pre_op_count, optimized_count

