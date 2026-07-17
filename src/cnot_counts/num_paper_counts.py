from src.optimize.paper_op import num_paper_gate_list, optimize_gate_list, count_cnots
from src.pauli_string_formation.mapping import pseudo_alphabetical_qubit_operator
from openfermion import QubitOperator

def num_paper_cnot_counts(op: QubitOperator):
    """
    Compute CNOT counts using the paper-style optimizer.

    Steps:
        1. Sort Pauli terms using pseudo-alphabetical ordering
        2. Convert to gate list (basis change + CNOT ladder + RZ + undo)
        3. Run custom optimizer (cancellation + merging + 3-CNOT rule)
        4. Count CNOTs before and after optimization

    Returns
    -------
    raw_cnot_count : int
        CNOT count before optimization
    opt_cnot_count : int
        CNOT count after optimization
    """
    ordered_terms = pseudo_alphabetical_qubit_operator(op)

    raw_gates = num_paper_gate_list(ordered_terms) # the actual list of gates before paper optimization
    opt_gates = optimize_gate_list(raw_gates)

    pre_op_count = count_cnots(raw_gates)
    optimized_count = count_cnots(opt_gates)

    return pre_op_count, optimized_count