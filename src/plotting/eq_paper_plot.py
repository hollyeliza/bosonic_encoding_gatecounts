import matplotlib.pyplot as plt

from pauli_string_formation.bosonic_disp import bosonic_disp_operator_matrix
from your_module_name import (
    naive_cnot_upper_bound_for_pair,
    diagonal_naive_cnot_upper_bound,
)


def run_counts(d_values):
    """
    Compute naive CNOT counts from the combinatorial equation
    for each encoding as a function of cutoff d.
    """
    encodings = ["sb", "gray", "unary"]
    results = {enc: [] for enc in encodings}

    for d in d_values:
        print(f"Running d={d}")
        qmat = bosonic_disp_operator_matrix(d)

        for enc in encodings:
            total = 0

            # Loop over Hermitian pairs (avoid double counting)
            for l in range(d):
                for lp in range(l, d):
                    coeff = qmat[l, lp]

                    if abs(coeff) < 1e-12:
                        continue

                    if l == lp:
                        total += diagonal_naive_cnot_upper_bound(l, d, enc)
                    else:
                        total += naive_cnot_upper_bound_for_pair(l, lp, d, enc)

            results[enc].append(total)
            print(f"  {enc:5s} | naive CX={total:6d}")

    return results


def plot_equation_results(d_values, data, title, filename):
    """
    Plot naive CNOT counts from the combinatorial equation.
    """
    plt.figure(figsize=(7, 4.5))

    styles = {
        "sb":   {"color": "#1c97ef", "linestyle": "-",  "marker": "o"},
        "gray": {"color": "#d7e727", "linestyle": "--", "marker": "s"},
        "unary":{"color": "#e66006", "linestyle": "-.", "marker": "^"},
    }

    labels = {
        "sb": "Std. binary",
        "gray": "Gray",
        "unary": "Unary"
    }

    for enc in ["sb", "gray", "unary"]:
        plt.plot(
            d_values,
            data[enc],
            label=labels[enc],
            **styles[enc]
        )

    plt.xlabel("Bosonic cutoff $d$")
    plt.ylabel("CNOT count")
    plt.title(title)

    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.savefig(filename, dpi=200)
    plt.show()


if __name__ == "__main__":
    d_values = list(range(2, 32))

    results = run_counts(d_values)

    plot_equation_results(
        d_values,
        results,
        title="Naive CNOT counts (combinatorial upper bound)",
        filename="results/plots/equation_cnot_counts.png",
    )