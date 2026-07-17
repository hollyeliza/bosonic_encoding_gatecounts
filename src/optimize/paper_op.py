from dataclasses import dataclass
from typing import Optional

# Attempt at construction of the optimizer used in the Sawaya et al. paper to obtain
# the CNOT counts.

# Start with op in pesudo alphabetical order
# construct CNOT ladder in this order
# make a class for the gate where attributes are: inverse gate, commutation properties with other gates, whether it is a rotation gate
# actual optimization 1:
# after contructing ladder, commute each gate as far forward and as far backwards as possible to see if you can cancel with inverse or merge
# actual optimization 2:
# pattern searching that cuts a 3 cnot gates into 2

@dataclass(frozen=True)
class Gate:
    """
    Minimal gate record for the SI-style circuit optimizer.

    name:
        "CNOT", "H", "B", "RZ"
        where:
          - H is a Hadamard
          - B is any self-inverse single-qubit basis-change gate
            such as the iR_{x+y}(pi)-type basis change mentioned in the SI
          - RZ is a Z rotation with angle theta
    qubits:
        tuple of qubit indices. For example:
          - H on qubit 3     -> (3,)
          - CNOT 1 -> 4      -> (1, 4)
    angle:
        only used for rotation gates such as RZ
    """
    name: str
    qubits: tuple[int, ...]
    angle: float = 0.0

    def acts_on(self) -> set[int]:
        return set(self.qubits)

    def is_rotation(self) -> bool:
        return self.name == "RZ"

    def is_self_inverse(self) -> bool:
        return self.name in {"CNOT", "H", "B"}

    def inverse(self) -> "Gate":
        if self.name in {"CNOT", "H", "B"}:
            return self
        if self.name == "RZ":
            return Gate("RZ", self.qubits, -self.angle)
        raise ValueError(f"No inverse rule for gate {self}")

    def same_support_and_type(self, other: "Gate") -> bool:
        return self.name == other.name and self.qubits == other.qubits


def count_cnots(gates: list[Gate]) -> int:
    """Count CNOT gates in a gate list."""
    return sum(g.name == "CNOT" for g in gates)


def remove_zero_rotations(gates: list[Gate], tol: float = 1e-12) -> list[Gate]:
    """Drop RZ gates with numerically zero angle."""
    out = []
    for g in gates:
        if g.name == "RZ" and abs(g.angle) < tol:
            continue
        out.append(g)
    return out


def commute_allowed(left: Gate, right: Gate) -> bool:
    """
    Conservative commutation rule.

    We only allow commuting when the supports are disjoint.
    This is safe and already captures much of the SI logic.

    You can later extend this with more specific rules if needed.
    """
    return left.acts_on().isdisjoint(right.acts_on())


def try_cancel_pair(a: Gate, b: Gate, tol: float = 1e-12) -> bool:
    """
    Return True if gates a and b cancel.
    """
    if a.name != b.name or a.qubits != b.qubits:
        return False

    if a.is_self_inverse():
        return True

    if a.is_rotation() and abs(a.angle + b.angle) < tol:
        return True

    return False


def try_merge_pair(a: Gate, b: Gate, tol: float = 1e-12) -> Optional[Gate]:
    """
    Merge adjacent rotations of the same type on the same qubit(s).
    """
    if not (a.is_rotation() and b.is_rotation()):
        return None
    if not a.same_support_and_type(b):
        return None

    merged_angle = a.angle + b.angle
    if abs(merged_angle) < tol:
        return Gate("RZ", a.qubits, 0.0)
    return Gate("RZ", a.qubits, merged_angle)


def bubble_once(records: list[Gate], start: int, direction: int) -> int:
    """
    Move one gate left (-1) or right (+1) through disjoint gates only.
    Returns the new index.
    """
    i = start
    while 0 <= i + direction < len(records):
        j = i + direction
        if commute_allowed(records[i], records[j]):
            records[i], records[j] = records[j], records[i]
            i = j
        else:
            break
    return i


def expose_cancellations_and_merges(records: list[Gate], tol: float = 1e-12) -> bool:
    """
    Simpler adjacent-only optimizer pass.

    It only checks gates that are already next to each other:
      - identical self-inverse gates cancel
      - opposite RZ rotations cancel
      - neighbouring RZ rotations merge
    """
    changed = False
    i = 0

    while i < len(records) - 1:
        left = records[i]
        right = records[i + 1]

        if try_cancel_pair(left, right, tol=tol):
            del records[i:i + 2]
            changed = True
            i = max(i - 1, 0)
            continue

        merged = try_merge_pair(left, right, tol=tol)
        if merged is not None:
            del records[i:i + 2]
            if abs(merged.angle) >= tol:
                records.insert(i, merged)
            changed = True
            i = max(i - 1, 0)
            continue

        i += 1

    return changed

def apply_three_cnot_rule(records: list[Gate]) -> bool:
    """
    Apply the explicit SI rewrite:

        CNOT(i0,i1) * CNOT(i1,i2) * CNOT(i0,i1)
            ->
        CNOT(i0,i2) * CNOT(i1,i2)

    Returns True if any replacement occurred.
    """
    changed = False
    i = 0

    while i <= len(records) - 3:
        a, b, c = records[i], records[i + 1], records[i + 2]

        if a.name == b.name == c.name == "CNOT":
            i0, i1 = a.qubits
            j0, j1 = b.qubits
            k0, k1 = c.qubits

            if (i0, i1) == (k0, k1) and i1 == j0:
                # pattern found: CNOT(i0,i1), CNOT(i1,i2), CNOT(i0,i1)
                i2 = j1
                replacement = [
                    Gate("CNOT", (i0, i2)),
                    Gate("CNOT", (i1, i2)),
                ]
                records[i:i + 3] = replacement
                changed = True
                i = max(i - 2, 0)
                continue

        i += 1

    return changed


def optimize_gate_list(gates, max_passes=50, tol=1e-12):
    records = list(gates)

    for pass_num in range(max_passes):
        before = list(records)
        before_len = len(records)

        changed = False
        changed |= expose_cancellations_and_merges(records, tol=tol)
        changed |= apply_three_cnot_rule(records)
        changed |= expose_cancellations_and_merges(records, tol=tol)

        records = remove_zero_rotations(records, tol=tol)

        print(f"pass {pass_num}: {before_len} -> {len(records)}", flush=True)

        if records == before:
            return records

        if len(records) > before_len:
            raise RuntimeError("Optimizer increased gate count; likely rewrite oscillation.")

    raise RuntimeError("Optimizer did not converge.")

def basis_change_gates_for_pauli(term) -> list[Gate]:
    """
    Return the forward basis-change gates needed to turn a Pauli string
    into a Z-string before the parity ladder.

    With the current abstract gate model:
      - X and Y are both represented by a generic self-inverse basis-change gate B
      - Z needs no basis change

    Important:
        This is good enough for CNOT-count comparisons under the paper-style
        optimizer, but it does NOT distinguish X-basis changes from Y-basis
        changes. So it is a simplified optimizer model, not a fully faithful
        single-qubit synthesis model.
    """
    gates: list[Gate] = []
    for q, p in term:
        if p in {"X", "Y"}:
            gates.append(Gate("B", (q,)))
        elif p == "Z":
            pass
        else:
            raise ValueError(f"Unexpected Pauli: {p}")
    return gates


def undo_basis_change_gates_for_pauli(term) -> list[Gate]:
    """
    Return the inverse basis-change gates.

    Since B is modeled as self-inverse, undoing is just the reverse list of B gates.
    """
    gates: list[Gate] = []
    for q, p in reversed(term):
        if p in {"X", "Y"}:
            gates.append(Gate("B", (q,)))
        elif p == "Z":
            pass
        else:
            raise ValueError(f"Unexpected Pauli: {p}")
    return gates


def pauli_term_to_gate_list(term, theta: float) -> list[Gate]:
    """
    Convert a single Pauli-string evolution exp(-i theta P) into the abstract
    gate list used by the paper-style optimizer.

    Parameters
    ----------
    term
        A Pauli term in OpenFermion/QubitOperator style, e.g.
        ((0, 'X'), (1, 'Y'), (3, 'Z'))
    theta
        Rotation parameter for exp(-i theta P)

    Returns
    -------
    list[Gate]
        Gate sequence:
            basis changes
            CNOT ladder
            RZ on target
            inverse CNOT ladder
            undo basis changes

    Notes
    -----
    - For a p-qubit Pauli string, this produces 2(p-1) CNOTs when p >= 2.
    - Identity terms produce no gates.
    """
    if len(term) == 0:
        return []

    gates: list[Gate] = []

    # 1. Basis changes so that P is turned into a Z-string
    gates.extend(basis_change_gates_for_pauli(term))

    qubits = [q for q, _ in term]

    # 2. Parity ladder + central RZ
    if len(qubits) == 1:
        gates.append(Gate("RZ", (qubits[0],), 2.0 * theta))
    else:
        target = qubits[-1]

        for control in qubits[:-1]:
            gates.append(Gate("CNOT", (control, target)))

        gates.append(Gate("RZ", (target,), 2.0 * theta))

        for control in reversed(qubits[:-1]):
            gates.append(Gate("CNOT", (control, target)))

    # 3. Undo basis changes
    gates.extend(undo_basis_change_gates_for_pauli(term))

    return gates


def num_paper_gate_list(op) -> list[Gate]:
    """
    Convert a QubitOperator into one flat gate list by exponentiating
    each Pauli term separately in the chosen order.

    Parameters
    ----------
    op
        QubitOperator in pseudoalphabetical order.

    Returns
    -------
    list[Gate]
        Concatenated gate list for all Pauli-string evolutions.

    Important
    ---------
    This assumes coefficients are real to numerical tolerance.
    If a coefficient has a significant imaginary part, an error is raised.
    """

    gates: list[Gate] = []

    for term, coeff in op.terms.items():
        coeff_c = complex(coeff)

        if abs(coeff_c.imag) > 1e-12:
            raise ValueError(f"Coefficient for term {term} has imaginary part: {coeff}")

        theta = float(coeff_c.real)
        gates.extend(pauli_term_to_gate_list(term, theta))

    return gates
