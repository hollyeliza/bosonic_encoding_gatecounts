# Bosonic Encoding Gate Counts

This project studies the CNOT cost of implementing a truncated bosonic
position-quadrature operator on qubits using different encodings:

- standard binary (`sb`)
- Gray code (`gray`)
- unary (`unary`)

The goal is to reproduce and explore the gate-count scaling trends discussed in
Sawaya et al., npj Quantum Information (2020).

## What This Does

- Builds truncated position and momentum quadrature matrices.
- Maps bosonic Fock states to qubits using standard binary, Gray, and unary
  encodings.
- Expands the encoded operators into Pauli strings with OpenFermion.
- Synthesizes Qiskit circuits from the Pauli-string operators.
- Counts CNOT gates before and after Qiskit transpiler optimization.
- Saves JSON results that can be inspected in the notebooks under `results/`.
- Includes a small validation notebook comparing encoded displacement circuits
  with `bosonic-qiskit` Wigner plots.

## Project Layout

```text
src/
  pauli_string_formation/
    encodings_b.py        # Fock-state to qubit encoding rules
    mapping.py            # Bosonic matrix elements to Pauli strings
    ps_quadratures.py     # Position and momentum quadrature matrices

  cnot_counts/
    eq_paper_counts.py    # Equation-style CNOT-count estimates
    num_paper_counts.py   # Numerical paper-style counting
    qiskit_counts.py      # Qiskit circuit CNOT counts

  optimize/
    paper_op.py           # Paper-inspired Pauli-string ordering/optimization
    qiskit_comp_op.py     # Qiskit circuit synthesis and transpilation helpers

  plotting/
    eq_paper_plot.py      # Generate equation-count result files
    num_paper_plot.py     # Generate numerical paper-count result files
    qiskit_plot.py        # Generate Qiskit-count result files

  bosonic_validation/
    bosonic_test_op.py    # Helpers for encoded displacement validation
    displacement.ipynb    # Wigner-plot validation notebook

results/
  *.json                  # Saved CNOT-count data
  *.ipynb                 # Plotting and inspection notebooks

tests/
  test_bosonic_disp.py
  test_encodings_b.py
  test_mapping.py
```

The `bosonic-qiskit/` directory is vendored for validation notebooks and Wigner
plotting experiments.

## Setup

From the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If you use the validation notebook, make sure the vendored `bosonic-qiskit`
package is available on the notebook path. The notebook handles this when run
from inside the repository.

## Generate Qiskit CNOT Counts

The Qiskit CNOT-count workflow starts in
`src/plotting/qiskit_plot.py`. This script builds the encoded displacement
operator, synthesizes the corresponding Qiskit circuit, counts CNOT gates before
and after Qiskit optimization, and saves the counts as JSON.

Before running it, edit the parameters in `src/plotting/qiskit_plot.py`,
especially:

- `d_values`
- `time`
- `trotter_steps`
- the output JSON filename

Then run this from the project root:

```bash
source .venv/bin/activate
python -m src.plotting.qiskit_plot
```

The script writes a JSON file to `results/`. The JSON stores, for each
encoding, the cutoff values and the CNOT counts before and after Qiskit
optimization.

To inspect the generated counts:

1. Open `results/qiskit_cnot.ipynb`.
2. Set `DATA_FILE` to the JSON filename you just generated, for example:

   ```python
   DATA_FILE = "qiskit_cnot_vs_d_150_16_07_trotter10.json"
   ```

3. Run the notebook cells to plot the raw and optimized CNOT counts.
4. Use the final printed values, or the values at a chosen cutoff such as
   `target_d = 128`, as the displacement CNOT counts for the Sawaya-style
   analysis.

For a quick check with an existing result file, open `results/qiskit_cnot.ipynb`
and use one of the JSON files already in `results/`.

## Sawaya-Style Analysis

The Sawaya notebooks use displacement CNOT counts to estimate the total
entangling-gate cost of a larger Hamiltonian simulation.

- `results/sawaya_ella.ipynb` is the original exploratory notebook. It contains
  hand-entered and WebPlotDigitizer data from Sawaya-style plots, plus several
  draft calculations for Hamiltonian gate-count estimates.
- `results/sawaya_ella_codex.ipynb` is the cleaned-up copy. It keeps the same
  idea but is organized top-to-bottom:
  - digitized Sawaya displacement CNOT counts
  - a plot recreating the Sawaya displacement-count curves
  - one parameter cell for the Hamiltonian estimate
  - helper functions for Trotter-step and entangling-gate estimates
  - separate comparisons using Sawaya counts, Qiskit 1-step displacement counts,
    and Qiskit 10-step displacement counts

The intended workflow is:

1. Generate Qiskit CNOT-count JSON with `python -m src.plotting.qiskit_plot`.
2. Open `results/qiskit_cnot.ipynb` and load that JSON with `DATA_FILE`.
3. Read off the optimized CNOT count for each encoding at the cutoff you want.
4. Put those displacement counts into `results/sawaya_ella_codex.ipynb`.
5. Run the Sawaya notebook to compare all-qubit encodings with the qubit-boson
   baseline.

## Useful Notebooks

- `results/qiskit_cnot.ipynb`: plots Qiskit CNOT counts from a saved JSON file.
- `results/sawaya_ella.ipynb`: original exploratory Sawaya comparison notebook.
- `results/sawaya_ella_codex.ipynb`: cleaned Sawaya comparison notebook for
  analyzing displacement CNOT counts inside the larger Hamiltonian estimate.

