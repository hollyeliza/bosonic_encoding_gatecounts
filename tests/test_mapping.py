import pytest
from openfermion import QubitOperator

from mapping import (
    one_qubit_map,
    matrix_element_to_qubit_operator,
    pauli_length,
    naive_cnot_count_from_qubit_operator,
)


def test_one_qubit_map_01():
    assert one_qubit_map(0, 1) == [("X", 0.5), ("Y", 0.5j)]


def test_one_qubit_map_10():
    assert one_qubit_map(1, 0) == [("X", 0.5), ("Y", -0.5j)]


def test_one_qubit_map_00():
    assert one_qubit_map(0, 0) == [("I", 0.5), ("Z", 0.5)]


def test_one_qubit_map_11():
    assert one_qubit_map(1, 1) == [("I", 0.5), ("Z", -0.5)]


def test_one_qubit_map_invalid():
    with pytest.raises(ValueError):
        one_qubit_map(2, 0)


def test_matrix_element_to_qubit_operator_unary_5_to_6():
    op = matrix_element_to_qubit_operator(5, 6, 1.0, 12, "unary")

    expected = (
        QubitOperator("X5 X6", 0.25)
        + QubitOperator("X5 Y6", 0.25j)
        + QubitOperator("Y5 X6", -0.25j)
        + QubitOperator("Y5 Y6", 0.25)
    )

    assert op == expected


def test_matrix_element_to_qubit_operator_unary_diagonal():
    op = matrix_element_to_qubit_operator(3, 3, 1.0, 8, "unary")

    expected = (
        QubitOperator((), 0.25)          # I on both qubits
        - QubitOperator("Z3", 0.25)
        + QubitOperator("Z3 Z3", 0.0)    # just to show intent; not needed
    )

    # Better to test term-by-term for unary diagonal support
    assert op.terms[((3, "Z"),)] == -0.5
    assert op.terms[()] == 0.5


def test_matrix_element_to_qubit_operator_gray_adjacent_levels_has_single_qubit_support(): # Failed
    op = matrix_element_to_qubit_operator(5, 6, 1.0, 16, "gray")

    # 5 -> 6 in Gray should differ by one bit only, so every non-identity term
    # should act on exactly one qubit.
    for term in op.terms:
        assert len(term) <= 1


def test_pauli_length():
    assert pauli_length(()) == 0
    assert pauli_length(((0, "X"),)) == 1
    assert pauli_length(((0, "X"), (2, "Y"))) == 2


def test_naive_cnot_count_single_two_qubit_term():
    op = QubitOperator("X0 Y1", 1.0)
    assert naive_cnot_count_from_qubit_operator(op) == 2


def test_naive_cnot_count_two_terms():
    op = QubitOperator("X0 Y1", 1.0) + QubitOperator("Z0 X1 Y2", 1.0)
    # lengths 2 and 3 -> 2*(2-1) + 2*(3-1) = 2 + 4 = 6
    assert naive_cnot_count_from_qubit_operator(op) == 6


def test_naive_cnot_count_ignores_single_qubit_and_identity_terms():
    op = QubitOperator((), 1.0) + QubitOperator("X0", 1.0) + QubitOperator("X0 Y1", 1.0)
    assert naive_cnot_count_from_qubit_operator(op) == 2