import matplotlib.pyplot as plt

from cnot_counts.qiskit_counts import run_counts 

def plot_qiskit_results(d_values, compiled, title, filename):
    """
    Plot compiled (Qiskit) CNOT counts vs cutoff d for all encodings.
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
            compiled[enc],
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


# ------------------ RUN ------------------

if __name__ == "__main__":
    d_values = list(range(2, 32))

    compiled = run_counts(d_values)

    plot_qiskit_results(
        d_values,
        compiled,
        title="Bosonic position operator $q$: compiled CNOT counts (Qiskit)",
        filename="results/plots/qiskit_compiled_cnot_counts.png",
    )