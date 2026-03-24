import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit

from bosonic_q import bosonic_q_matrix
from mapping import (
    matrix_to_qubit_operator,
    naive_cnot_count_from_qubit_operator,
    sorted_terms_pseudo_alphabetical,
)
from encodings import n_qubits
from synthesize import append_pauli_evolution, cnot_count_of_compiled_circuit

def build_full_circuit_from_operator(op, nq: int, theta_scale: float = 1.0):
    """
    Build a circuit by concatenating exponentials of Pauli terms
    in pseudo-alphabetical order. For gate counting only.
    """
    qc = QuantumCircuit(nq)
    for term, coeff in sorted_terms_pseudo_alphabetical(op):
        if abs(coeff) < 1e-12:
            continue
        theta = float(np.real(coeff)) * theta_scale
        # For q-hat all coefficients should be real after Hermitian summation.
        append_pauli_evolution(qc, term, theta)
    return qc


def run_counts(d_values):
    encodings = ["sb", "gray", "unary"]

    naive = {enc: [] for enc in encodings}
    compiled = {enc: [] for enc in encodings}

    for d in d_values:
        print(f"Running d={d}")
        qmat = bosonic_q_matrix(d)

        for enc in encodings:
            op = matrix_to_qubit_operator(qmat, d=d, encoding=enc)
            naive_count = naive_cnot_count_from_qubit_operator(op)

            nq = n_qubits(d, enc)
            qc = build_full_circuit_from_operator(op, nq=nq, theta_scale=1.0)
            compiled_count, tqc = cnot_count_of_compiled_circuit(qc, optimization_level=3)

            naive[enc].append(naive_count)
            compiled[enc].append(compiled_count)

            print(
                f"  {enc:5s} | nq={nq:2d} | #terms={len(op.terms):4d} | "
                f"naive CX={naive_count:5d} | compiled CX={compiled_count:5d}"
            )
    return naive, compiled


def plot_results(d_values, data, title, filename):
    plt.figure(figsize=(7, 4.5))
    for enc, label in [("sb", "Std. binary"), ("gray", "Gray"), ("unary", "Unary")]:
        plt.plot(d_values, data[enc], marker="o", label=label)
    plt.xlabel("d")
    plt.ylabel("CNOT count")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.show()


if __name__ == "__main__":
    d_values = list(range(2, 17))
    naive, compiled = run_counts(d_values)

    plot_results(
        d_values,
        naive,
        title="Bosonic position operator q: naive CNOT counts",
        filename="q_naive_cnot_counts.png",
    )
    plot_results(
        d_values,
        compiled,
        title="Bosonic position operator q: compiled CNOT counts",
        filename="q_compiled_cnot_counts.png",
    )