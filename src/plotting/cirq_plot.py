"""
Plot Cirq-optimized CNOT cost vs cutoff d for different encodings.

This script assumes your codebase already provides:
  - bosonic_disp_operator_matrix(d)
  - matrix_to_qubit_operator(mat, d, encoding, tol=...)
  - cirq_cnot_count_after_optimization(op, time=..., term_order=...)

Edit the import block below to match your repo layout if needed.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt

from src.pauli_string_formation.mapping import matrix_to_qubit_operator, pseudo_alphabetical_qubit_operator, bosonic_disp_operator_matrix
from src.cnot_counts.cirq_counts import cirq_cnot_count_before_and_after_optimization

from collections.abc import Iterable
from pathlib import Path
import matplotlib.pyplot as plt


def cirq_counts_for_encoding(
    d_values: Iterable[int],
    encoding: str,
) -> tuple[list[int], list[int], list[int]]:

    """
    For each d (we want to calculate the number of CNOTs for different bosonic truncations):
      1. build q operator matrix
      2. map to qubit operator under chosen encoding
      3. count raw and Cirq-optimized CNOTs

    Returns
    -------
    ds, raw_counts, opt_counts
    """
    ds: list[int] = []
    pre_op_count: list[int] = []
    optimized_count: list[int] = []

    for d in d_values:
        print(f"[{encoding}] building operator for d={d}")
        mat = bosonic_disp_operator_matrix(d)
        op = matrix_to_qubit_operator(mat, d, encoding)
        op_pseudo = pseudo_alphabetical_qubit_operator(op)

        raw_cx, opt_cx = cirq_cnot_count_before_and_after_optimization(op_pseudo) 

        print(f"[{encoding}] d={d}: raw={raw_cx}, optimized={opt_cx}")

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
    Optionally overlay raw counts as dashed lines.
    """
    plt.figure(figsize=(8, 5))

    for encoding, (ds, pre_op_count, optimized_count) in results.items():
        plt.plot(ds, optimized_count, marker="o", label=f"{encoding} (optimized)")

    plt.xlabel("Cutoff d")
    plt.ylabel("CNOT count")
    plt.title("Cirq-optimized CNOT cost vs cutoff d")
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
        results[encoding] = cirq_counts_for_encoding(
            d_values=d_values,
            encoding=encoding
        )

    save_plot(
        results=results,
        output_path=Path("cirq_cnot_vs_d.png"),
        show_raw=True,
        logy=False,
    )

if __name__ == "__main__":
    main()