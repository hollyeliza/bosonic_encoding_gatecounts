# Bosonic Encoding Gate Counts

Reproduction of CNOT gate count scaling for the bosonic position operator
across different encodings (standard binary, Gray, unary), following:

Sawaya et al., npj Quantum Information (2020)

The unoptimized CNOT counts are looked at alongside different optimization methods
(the way the paper described and also qiskit and cirq) to try and recreate the CNOT
counts obtained in the paper for the encodings.

Refer to the Wiki.

## What this does

- Constructs truncated bosonic displacement operator
- Maps it to qubit operators using different encodings
- Expands into Pauli strings using OpenFermion
- Synthesizes circuits and counts CNOT gates (naive equation way and using compilation)
- 'Tries to' reproduce scaling trends in the paper


## File structure

Note: Since all changed - need to update

src/
- pauli_string_formation.py # maps bosonic displacement operator to sum of pauli strings for different encodings (unary, Gray, binary)
- cnot_counts.py            # counts no. of entangling (cnot) gates to implement bosonic displacement operator on all-qubit quantum computer
- optimize.py        # optimization procedures to cut down the number of cnot gates (qiskit, cirq, paper rationale)
- plotting.py         # reproduces plots for CNOT count vs phonon occupation number

tests/
- test_mapping.py       # unit tests for mapping functions


## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
