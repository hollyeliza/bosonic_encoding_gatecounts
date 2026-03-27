# This is currently specific to qiskit
# Make specific - make sure you sort alphabetical and prep according to SI before passing
# into the optimizer

import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit

from bosonic_disp import bosonic_disp_operator_matrix
from mapping import (
    matrix_to_qubit_operator,
    naive_cnot_count_from_qubit_operator,
    sorted_terms_pseudo_alphabetical,
)
from encodings_b import n_qubits
from qiskit_comp_op import append_pauli_evolution, cnot_count_of_compiled_circuit

def build_full_circuit_from_operator(op, nq: int, theta_scale: float = 1.0):
    """
    Construct a circuit for a Pauli-operator decomposition by appending
    exp(-i θ_k P_k) for each Pauli term P_k in a fixed order.

    Here `op` is a QubitOperator of the form
        op = sum_k coeff_k * P_k.

    For each term:
        θ_k = coeff_k * theta_scale

    The resulting circuit is intended for gate-count estimation rather than
    high-accuracy simulation, so only the real part of each coefficient is used.

    Warning: Not sure about QubitOperator implementation.
    """
    qc = QuantumCircuit(nq)
    for term, coeff in sorted_terms_pseudo_alphabetical(op):
        if abs(coeff) < 1e-12:
            continue
        theta = float(np.real(coeff)) * theta_scale # This sets the evolution angle for that Pauli term.
        # For q-hat all coefficients should be real after Hermitian summation.
        append_pauli_evolution(qc, term, theta)
    return qc


def run_counts(d_values):
    """
    For each cutoff d in d_values and for each encoding ("sb", "gray", "unary"),
    compute:

    1. the encoded Pauli representation of the bosonic position operator q,
    2. the naive CNOT count using the rule 2(p-1) for each Pauli string,
    3. the compiled CNOT count obtained by building and transpiling a circuit.

    Returns
    -------
    naive : dict[str, list[int]]
        Naive CNOT counts for each encoding as a function of d.
    compiled : dict[str, list[int]]
        Transpiled CNOT counts for each encoding as a function of d.
    """
    encodings = ["sb", "gray", "unary"]

    naive = {enc: [] for enc in encodings} # using 2p - 2 from paper
    compiled = {enc: [] for enc in encodings} # using qiskit

    for d in d_values:
        print(f"Running d={d}")
        qmat = bosonic_disp_operator_matrix(d)

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
    """"
    Plot CNOT count versus bosonic cutoff d for the three encodings.
    """
    plt.figure(figsize=(7, 4.5))
    styles = {
    "sb":   {"color": "#1c97ef", "linestyle": "-",  "marker": "o"},
    "gray": {"color": "#d7e727", "linestyle": "--", "marker": "s"},
    "unary":{"color": "#e66006", "linestyle": "-.", "marker": "^"},
    }

    for enc, label in [("sb", "Std. binary"), ("gray", "Gray"), ("unary", "Unary")]:
        plt.plot(
            d_values,
            data[enc],
            label=label,
            **styles[enc]
        )
    plt.xlabel("d")
    plt.ylabel("CNOT count")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.show()


def plot_combined_results(d_values, naive, compiled, title, filename):
    """
    Plot naive and compiled CNOT counts on the same axes.
    Same color = same encoding
    Different linestyle = naive vs compiled
    """
    plt.figure(figsize=(7, 4.5))

    styles = {
        "sb":   {"color": "#1c97ef"},
        "gray": {"color": "#d7e727"},
        "unary":{"color": "#e66006"},
    }

    for enc, label in [("sb", "Std. binary"), ("gray", "Gray"), ("unary", "Unary")]:

        # naive → dashed
        plt.plot(
            d_values,
            naive[enc],
            linestyle="--",
            marker="o",
            label=f"{label} (naive)",
            **styles[enc]
        )

        # compiled → solid
        plt.plot(
            d_values,
            compiled[enc],
            linestyle="-",
            marker="o",
            label=f"{label} (compiled)",
            **styles[enc]
        )

    plt.xlabel("d")
    plt.ylabel("CNOT count")
    plt.title(title)

    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.savefig(filename, dpi=200)
    plt.show()


if __name__ == "__main__":
    d_values = list(range(2, 32))
    naive, compiled = run_counts(d_values)

    # Naive (logic and equation in the paper - 2(p-1) for each Pauli String)
    # plot_results(
    #     d_values,
    #     naive,
    #     title="Bosonic position operator q: naive CNOT counts",
    #     filename="q_naive_cnot_counts.png",
    # )

    # Compiled (qiskit)
    plot_results(
        d_values,
        compiled,
        title="Bosonic position operator q: compiled CNOT counts",
        filename="q_compiled_cnot_counts.png",
    )

    # Combined comparison
    # plot_combined_results(
    #     d_values,
    #     naive,
    #     compiled,
    #     title="Bosonic position operator q: naive vs compiled CNOT counts",
    #     filename="q_combined_cnot_counts.png",
    # )

    # get compiled counts