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

try:
    from optimize.cirq_comp_op import cirq_cnot_count_after_optimization
except Exception:
    # Fallback if you put the function in a different module
    from optimize.cirq_comp_op import (
        qubit_operator_to_trotter_circuit,
        optimize_cirq_circuit,
        count_cnots,
    )
    from openfermion import QubitOperator
    import cirq

    def cirq_cnot_count_after_optimization(
        op: QubitOperator,
        time: float = 1.0,
        term_order: str = "default",
    ) -> tuple[int, int, cirq.Circuit, cirq.Circuit]:
        raw_circuit, _ = qubit_operator_to_trotter_circuit(
            op,
            time=time,
            term_order=term_order,
        )
        optimized_circuit = optimize_cirq_circuit(raw_circuit)
        raw_count = count_cnots(raw_circuit)
        optimized_count = count_cnots(optimized_circuit)
        return raw_count, optimized_count, raw_circuit, optimized_circuit

from pauli_string_formation.mapping import matrix_to_qubit_operator, bosonic_disp_operator_matrix

def parse_d_values(d_arg: str) -> list[int]:
    """
    Parse either:
      - a comma-separated list, e.g. '2,3,4,5,6,7,8'
      - a range spec start:end, e.g. '2:16' -> [2, 3, ..., 16]
    """
    d_arg = d_arg.strip()
    if ":" in d_arg:
        start_str, end_str = d_arg.split(":", 1)
        start = int(start_str)
        end = int(end_str)
        if end < start:
            raise ValueError("For d range start:end, need end >= start.")
        return list(range(start, end + 1))
    return [int(x.strip()) for x in d_arg.split(",") if x.strip()]


def cirq_counts_for_encoding(
    d_values: Iterable[int],
    encoding: str,
    tol: float = 1e-14,
    time: float = 1.0,
    term_order: str = "pseudo_alpha",
) -> tuple[list[int], list[int], list[int]]:
    """
    For each d:
      1. build q operator matrix
      2. map to qubit operator under chosen encoding
      3. count raw and Cirq-optimized CNOTs

    Returns
    -------
    ds, raw_counts, opt_counts
    """
    ds: list[int] = []
    raw_counts: list[int] = []
    opt_counts: list[int] = []

    for d in d_values:
        print(f"[{encoding}] building operator for d={d}")
        mat = bosonic_disp_operator_matrix(d)
        op = matrix_to_qubit_operator(mat, d, encoding, tol=tol)

        raw_cx, opt_cx, _, _ = cirq_cnot_count_after_optimization(
            op,
            time=time,
            term_order=term_order,
        )

        print(f"[{encoding}] d={d}: raw={raw_cx}, optimized={opt_cx}")

        ds.append(d)
        raw_counts.append(raw_cx)
        opt_counts.append(opt_cx)

    return ds, raw_counts, opt_counts


def save_plot(
    results: dict[str, tuple[list[int], list[int], list[int]]],
    output_path: Path,
    show_raw: bool = False,
    logy: bool = False,
) -> None:
    """
    Plot optimized CNOT count vs d for each encoding.
    Optionally also overlay raw counts.
    """
    plt.figure(figsize=(8, 5))

    for encoding, (ds, raw_counts, opt_counts) in results.items():
        plt.plot(ds, opt_counts, marker="o", label=f"{encoding} (optimized)")
        if show_raw:
            plt.plot(ds, raw_counts, marker="x", linestyle="--", label=f"{encoding} (raw)")

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
    parser = argparse.ArgumentParser(
        description="Plot Cirq-optimized CNOT cost vs cutoff d for different encodings."
    )
    parser.add_argument(
        "--encodings",
        type=str,
        default="sb,gray,unary",
        help="Comma-separated encodings to compare, e.g. 'sb,gray,unary'.",
    )
    parser.add_argument(
        "--d-values",
        type=str,
        default="2:16",
        help="Either comma-separated values like '2,3,4,5' or range '2:16'.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-14,
        help="Tolerance passed to matrix_to_qubit_operator.",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=1.0,
        help="Evolution time passed to cirq_cnot_count_after_optimization.",
    )
    parser.add_argument(
        "--term-order",
        type=str,
        default="pseudo_alpha",
        choices=["default", "pseudo_alpha"],
        help="Pauli term ordering for circuit construction.",
    )
    parser.add_argument(
        "--show-raw",
        action="store_true",
        help="Overlay raw (pre-optimization) CNOT counts.",
    )
    parser.add_argument(
        "--logy",
        action="store_true",
        help="Use a log scale on the y-axis.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="cirq_cnot_vs_d.png",
        help="Output plot filename.",
    )

    args = parser.parse_args()

    encodings = [x.strip() for x in args.encodings.split(",") if x.strip()]
    d_values = parse_d_values(args.d_values)

    results: dict[str, tuple[list[int], list[int], list[int]]] = {}
    for encoding in encodings:
        results[encoding] = cirq_counts_for_encoding(
            d_values=d_values,
            encoding=encoding,
            tol=args.tol,
            time=args.time,
            term_order=args.term_order,
        )

    save_plot(
        results=results,
        output_path=Path(args.output),
        show_raw=args.show_raw,
        logy=args.logy,
    )


if __name__ == "__main__":
    main()
