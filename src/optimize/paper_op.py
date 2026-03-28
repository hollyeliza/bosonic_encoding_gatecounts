# code matching paper's logic:

# The circuit optimization was performed as follows. In the optimizer, every gate is represented 
# in a data structure that contains its attributes, e.g. its inverse gate, its commutation properties
# with other gates, and whether it is a rotation gate. For each gate, the optimizer looks for an 
# opportunity to cancel it with its inverse or merge it with another rotation gate by commuting it as 
# far forwards and backwards as possible. The optimizer also looks for a limited set of commonly 
# occurring patterns that allow for gate reductions. For the circuits used in this study, this 
# pattern-searching allows us to reduce some sets of three CNOT gates to two: any gate sequence 
# CNOT(i0,i1) × CNOT(i1 ,i2 ) × CNOT(i0 ,i1 ) is changed to the equivalent CNOT(i0 ,i2 ) × CNOT(i1 ,i2 ). 
# Merging and cancelling gates, as well as this pattern-searching, are performed for several passes until 
# the circuit converges. The choice of ordering for the Pauli terms affects the gate counts, as different 
# orderings of CNOT ladders affect the presence of particular pairs of eliminable gates. Finding the absolute 
# optimal ordering is a combinatorially hard problem and is not the focus of this work [BKM18]. However, to 
# test the quality of the default ordering, we took the encoded bosonic qˆ operator for the unary, gray, and 
# standard binary codes, testing >900 random orderings for each.
#
# We found that the default ordering was superior to every random ordering. This result is intuitive, as 
# orderings that match similar Pauli strings side-by-side will lead to more pairs of adjacent CNOT gates 
# and adjacent single-qubit gates.

# Convert Pauli term to CNOT ladder in order obtained after pseudo-alphabetical def
# Construct circuit
# Commut each term as far back and as far forward as possible to see if it can be cancelled with
# inverse or merged with another rotation gate.
# The CNOT pattern described is looked for
# This should be it

from dataclasses import dataclass
from typing import Optional


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
    For each gate, commute it backwards through disjoint gates as far as possible,
    then look for local cancellation / merging opportunities.

    Returns True if anything changed.
    """
    changed = False
    i = 0

    while i < len(records):
        # commute backwards as far as possible
        new_i = bubble_once(records, i, -1)
        if new_i != i:
            changed = True
            i = new_i

        local_changed = True
        while local_changed:
            local_changed = False

            # cancel with left neighbor
            if i > 0 and try_cancel_pair(records[i - 1], records[i], tol=tol):
                del records[i]
                del records[i - 1]
                changed = True
                local_changed = True
                i = max(i - 2, 0)
                continue

            # merge with left neighbor
            if i > 0:
                merged = try_merge_pair(records[i - 1], records[i], tol=tol)
                if merged is not None:
                    del records[i]
                    del records[i - 1]
                    if abs(merged.angle) >= tol:
                        records.insert(i - 1, merged)
                    changed = True
                    local_changed = True
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


def optimize_gate_list(
    gates: list[Gate],
    max_passes: int = 50,
    tol: float = 1e-12,
) -> list[Gate]:
    """
    SI-style repeated optimization until convergence.

    Each pass:
      1. commute gates to expose local simplifications
      2. cancel inverse pairs
      3. merge rotations
      4. apply the 3-CNOT -> 2-CNOT pattern
      5. repeat until no further change
    """
    records = list(gates)

    for _ in range(max_passes):
        changed = False

        if expose_cancellations_and_merges(records, tol=tol):
            changed = True

        if apply_three_cnot_rule(records):
            changed = True

        if expose_cancellations_and_merges(records, tol=tol):
            changed = True

        records = remove_zero_rotations(records, tol=tol)

        if not changed:
            break

    return records

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


def qubit_operator_to_gate_list(op, term_order=None, tol: float = 1e-12) -> list[Gate]:
    """
    Convert a QubitOperator into one flat gate list by exponentiating
    each Pauli term separately in the chosen order.

    Parameters
    ----------
    op
        QubitOperator-like object with .terms dictionary
    term_order
        Optional iterable of (term, coeff) pairs. If omitted, uses op.terms.items().
        This lets you pass in your pseudo-alphabetically sorted ordering.
    tol
        Numerical tolerance for dropping tiny coefficients.

    Returns
    -------
    list[Gate]
        Concatenated gate list for all Pauli-string evolutions.

    Important
    ---------
    This assumes coefficients are real to numerical tolerance.
    If a coefficient has a significant imaginary part, an error is raised.
    """
    if term_order is None:
        term_order = op.terms.items()

    gates: list[Gate] = []

    for term, coeff in term_order:
        if abs(coeff) < tol:
            continue

        if len(term) == 0:
            # ignore identity term for gate counting
            continue

        coeff_c = complex(coeff)
        if abs(coeff_c.imag) > tol:
            raise ValueError(
                f"Non-negligible imaginary coefficient encountered for term {term}: {coeff}"
            )

        theta = float(coeff_c.real)
        gates.extend(pauli_term_to_gate_list(term, theta))

    return gates

