from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import json

from src.pauli_string_formation.mapping import (
    matrix_to_qubit_operator,
    pseudo_alphabetical_qubit_operator,
    bosonic_disp_operator_matrix,
)
from src.cnot_counts.qiskit_counts import qiskit_cnot_count_before_and_after_optimization


def qiskit_counts_for_encoding(
    d_values: Iterable[int],
    encoding: str,
) -> tuple[list[int], list[int], list[int]]:
    """
    For each d:
      1. build q operator matrix
      2. map to qubit operator under chosen encoding
      3. count raw and Qiskit-optimized CNOTs

    Returns
    -------
    ds, raw_counts, optimized_counts
    """
    ds: list[int] = []
    pre_op_count: list[int] = []
    optimized_count: list[int] = []

    for d in d_values:
        print(f"[{encoding}] building operator for d={d}")

        mat = bosonic_disp_operator_matrix(d)
        op = matrix_to_qubit_operator(mat, d, encoding)
        op_pseudo = pseudo_alphabetical_qubit_operator(op)

        raw_cx, opt_cx = qiskit_cnot_count_before_and_after_optimization(op_pseudo)

        print(f"[{encoding}] d={d}: pre-qiskit op cnot count={raw_cx}, post-qiskit op cnot count={opt_cx}")

        ds.append(d)
        pre_op_count.append(raw_cx)
        optimized_count.append(opt_cx)

    return ds, pre_op_count, optimized_count


def save_plot(
    results: dict[str, tuple[list[int], list[int], list[int]]],
    output_path: Path,
    show_raw: bool = False,
    logy: bool = False,
) -> None:
    """
    Plot optimized CNOT count vs cutoff d for each encoding.
    Optionally overlay raw counts.
    """
    plt.figure(figsize=(8, 5))

    for encoding, (ds, pre_op_count, optimized_count) in results.items():
        plt.plot(ds, optimized_count, marker="o", label=f"{encoding} (after optimization)")
        if show_raw:
            plt.plot(ds, pre_op_count, marker="o", linestyle="--", label=f"{encoding} (before optimization)")

    plt.xlabel("Cutoff d")
    plt.ylabel("CNOT count")
    plt.title("Qiskit pre and post optimized CNOT cost vs cutoff d")

    if logy:
        plt.yscale("log")

    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved plot to {output_path}")


def main() -> None:
    encodings = ["sb", "gray", "unary"]
    d_values = list(range(2, 17))

    results: dict[str, tuple[list[int], list[int], list[int]]] = {}

    for encoding in encodings:
        results[encoding] = qiskit_counts_for_encoding(
            d_values=d_values,
            encoding=encoding,
        )

    results_dir = Path(__file__).resolve().parents[1] / "results"
    results_dir.mkdir(exist_ok=True)

    # json for notebook
    json_path = results_dir / "qiskit_cnot_vs_d.json"
    with open(json_path, "w") as f:
        json.dump(results, f)

    # png
    plot_path = results_dir / "qiskit_cnot_vs_d.png"
    
    save_plot(
        results=results,
        output_path=plot_path,
        show_raw=True,
        logy=False,
    )

    save_plot(
        results=results,
        output_path=plot_path,
        show_raw=True,
        logy=False,
    )


if __name__ == "__main__":
    main()