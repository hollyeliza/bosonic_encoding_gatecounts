# Bosonic Encoding Gate Counts

Reproduction of CNOT gate count scaling for the bosonic position operator
across different encodings (standard binary, Gray, unary), following:

Sawaya et al., npj Quantum Information (2020)

## What this does

- Constructs truncated bosonic displacement operator
- Maps it to qubit operators using different encodings
- Expands into Pauli strings using OpenFermion
- Synthesizes circuits and counts CNOT gates (naive equation way and using compilation)
- 'Tries to' reproduce scaling trends in the paper

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt