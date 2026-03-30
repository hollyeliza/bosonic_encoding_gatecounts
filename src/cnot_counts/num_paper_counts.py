from optimize.paper_op import qubit_operator_to_gate_list, optimize_gate_list, count_cnots
from pauli_string_formation.mapping import sorted_terms_pseudo_alphabetical
from openfermion import QubitOperator

def paper_cnot_counts(op: QubitOperator):
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
    raw_gates : list[Gate]
        Unoptimized gate list
    opt_gates : list[Gate]
        Optimized gate list
    """
    ordered_terms = sorted_terms_pseudo_alphabetical(op)

    raw_gates = qubit_operator_to_gate_list(op, term_order=ordered_terms)
    opt_gates = optimize_gate_list(raw_gates)

    raw_cnot_count = count_cnots(raw_gates)
    opt_cnot_count = count_cnots(opt_gates)

    return raw_cnot_count, opt_cnot_count, raw_gates, opt_gates