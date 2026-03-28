import matplotlib.pyplot as plt

from cnot_counts.num_paper_counts import paper_cnot_counts


def plot_paper_results(d_values, raw, optimized, title, filename):
    """
    Plot paper-style CNOT counts (before vs after optimization).

    Same color = same encoding
    Different linestyle = raw vs optimized
    """
    plt.figure(figsize=(7, 4.5))

    colors = {
        "sb": "#1c97ef",
        "gray": "#d7e727",
        "unary": "#e66006",
    }

    labels = {
        "sb": "Std. binary",
        "gray": "Gray",
        "unary": "Unary",
    }

    for enc in ["sb", "gray", "unary"]:

        # raw → dashed
        plt.plot(
            d_values,
            raw[enc],
            linestyle="--",
            marker="o",
            color=colors[enc],
            label=f"{labels[enc]} (raw)"
        )

        # optimized → solid
        plt.plot(
            d_values,
            optimized[enc],
            linestyle="-",
            marker="o",
            color=colors[enc],
            label=f"{labels[enc]} (optimized)"
        )

    plt.xlabel("Bosonic cutoff $d$")
    plt.ylabel("CNOT count")
    plt.title(title)

    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.savefig(filename, dpi=200)
    plt.show()