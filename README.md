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
    qiskit_counts.py      # Qiskit circuit CNOT counts

  optimize/
    qiskit_comp_op.py     # Qiskit circuit synthesis and transpilation helpers

  plotting/
    qiskit_plot.py        # Generate Qiskit-count result files

  bosonic_validation/
    bosonic_test_op.py    # Helpers for encoded displacement validation
    displacement.ipynb    # Wigner-plot validation notebook
```

## Setup

From the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If you use the validation notebook, clone `bosonic-qiskit` into the project
root first:

```bash
git clone https://github.com/C2QA/bosonic-qiskit.git
```

The validation notebook adds `bosonic-qiskit/src` to the notebook path when run
from inside this repository.

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
