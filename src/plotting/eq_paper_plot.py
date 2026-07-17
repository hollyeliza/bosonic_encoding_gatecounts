import matplotlib.pyplot as plt
from collections.abc import Iterable
from pathlib import Path
import json


from pauli_string_formation.ps_quadratures import position_operator_matrix
from src.cnot_counts.eq_paper_counts import naive_cnot_upper_bound_for_transition


def equation_paper_counts_for_encoding(
    d_values: Iterable[int],
    encoding: str,
) -> tuple[list[int], list[int], list[int]]:
    """
    Same structure as qiskit_counts_for_encoding,
    but only naive counts exist.
    """
    ds: list[int] = []
    naive_counts: list[int] = []
    dummy_optimized: list[int] = []  # keep structure consistent

    for d in d_values:
        print(f"[{encoding}] computing equation counts for d={d}")

        qmat = position_operator_matrix(d)

        total = 0
        for l in range(d):
            for lp in range(l, d):
                coeff = qmat[l, lp]

                if abs(coeff) < 1e-12:
                    continue

                total += naive_cnot_upper_bound_for_transition(l, lp, d, encoding)

        print(f"[{encoding}] d={d}: naive CX={total}")

        ds.append(d)
        naive_counts.append(total)
        dummy_optimized.append(0)  # placeholder

    return ds, naive_counts, dummy_optimized


def save_plot(
    results: dict[str, tuple[list[int], list[int], list[int]]],
    output_path: Path,
    show_raw: bool = False,
    logy: bool = False,
) -> None:
    """
    Plot equation-based CNOT count vs cutoff d for each encoding.

    The equation count is the naive combinatorial upper bound from the paper.
    There is no optimized count here, so show_raw is kept only to match the
    plotting function structure used in the Qiskit plot file.
    """
    plt.figure(figsize=(8, 5))

    for encoding, (ds, naive_counts, _) in results.items():
        plt.plot(ds, naive_counts, marker="o", label=f"{encoding} (equation)")

    plt.xlabel("Cutoff d")
    plt.ylabel("CNOT count")
    plt.title("Naive CNOT cost vs cutoff d (combinatorial equation)")

    if logy:
        plt.yscale("log")

    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved plot to {output_path}")


def main() -> None:
    encodings = ["sb", "gray", "unary"]
    d_values = list(range(2, 33))

    results: dict[str, tuple[list[int], list[int], list[int]]] = {}

    for encoding in encodings:
        results[encoding] = equation_paper_counts_for_encoding(
            d_values=d_values,
            encoding=encoding,
        )

    results_dir = Path(__file__).resolve().parents[2] / "results"
    results_dir.mkdir(exist_ok=True)

    json_path = results_dir / "equation_cnot_vs_d_32_01_05.json"
    with open(json_path, "w") as f:
        json.dump(results, f)

    print(f"Saved JSON to {json_path}")

    plot_path = results_dir / "equation_cnot_vs_d_32_01_05.png"
    save_plot(
        results=results,
        output_path=plot_path,
        show_raw=False,
        logy=False,
    )


if __name__ == "__main__":
    main()