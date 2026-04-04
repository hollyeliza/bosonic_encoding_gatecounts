# Bosonic Encoding Gate Counts

Reproduction of CNOT gate count scaling for the bosonic position operator
across different encodings (standard binary, Gray, unary), following:

Sawaya et al., npj Quantum Information (2020)

The unoptimized CNOT counts are looked at alongside different optimization methods
(the way the paper described and also qiskit and cirq) to try and recreate the CNOT
counts obtained in the paper for the encodings.

## What this does

- Constructs truncated bosonic displacement operator
- Maps it to qubit operators using different encodings
- Expands into Pauli strings using OpenFermion
- Synthesizes circuits and counts CNOT gates (naive equation way and using compilation)
- 'Tries to' reproduce scaling trends in the paper


## File structure

Note: Since all changed - need to update

src/
- encodings_b.py        # bit encodings (unary, Gray, binary)
- mapping.py            # maps |l><l'| → Pauli strings
- bosonic_ops.py        # constructs bosonic operators (e.g. q)
- gate_count.py         # CNOT counting utilities

tests/
- test_mapping.py       # unit tests for mapping functions


## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
