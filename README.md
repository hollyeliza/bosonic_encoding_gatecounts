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

## To do

1. Equation in paper cnot counts and plot

2. Optimization in paper cnot counts and plot

3. I remove the double counting of the diagonal but there didn't seem to be a change in the counting

4. Write out optimisation process for cirq in wiki and the entire workflow to get the counts (will repeat same pattern for qiskit and papers)

- unsure if the crirq optimization has done anything as the count before and after seems to be the same ...

5. Cleanly write out the intuition of the CNOT ladder

6. Fully write out the wiki

7. Write tests

8. Pylint

9. Clean up the README (use it for the paper)

## File structure

Note: Since all changed - need to update

src/
- encodings_b.py        # bit encodings (unary, Gray, binary)
- mapping.py            # maps |l><l'| → Pauli strings
- bosonic_ops.py        # constructs bosonic operators (e.g. q)
- gate_count.py         # CNOT counting utilities

tests/
- test_mapping.py       # unit tests for mapping functions


## Key idea

Bosonic operators are sparse in the Fock basis but become non-local 
when mapped to qubits. This project quantifies the resulting gate overhead.


## Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
