import pytest

from pauli_string_formation.encodings_b import (
    n_qubits_sb,
    sb_bits,
    gray_int,
    gray_bits,
    unary_bits,
    bitmask_subset,
    bits_for_level,
    n_qubits,
)

def test_n_qubits_sb():
    assert n_qubits_sb(1) == 0
    assert n_qubits_sb(2) == 1
    assert n_qubits_sb(3) == 2
    assert n_qubits_sb(8) == 3
    assert n_qubits_sb(10) == 4
    assert n_qubits_sb(32) == 5


def test_sb_bits():
    assert sb_bits(0, 8) == [0, 0, 0]
    assert sb_bits(1, 8) == [1, 0, 0]
    assert sb_bits(8, 9) == [0, 0, 0, 1]
    assert sb_bits(3, 32) == [1, 1, 0, 0, 0] # [00011] in reality but read left to right - most important bit (representing 2^0) usually far right


def test_sb_bits_invalid_level():
    with pytest.raises(ValueError):
        sb_bits(5, 4)
        sb_bits(5, 5)


def test_gray_int():
    assert gray_int(0) == 0
    assert gray_int(1) == 1
    assert gray_int(2) == 3
    assert gray_int(3) == 2
    assert gray_int(4) == 6


def test_gray_bits():
    assert gray_bits(0, 8) == [0, 0, 0]
    assert gray_bits(1, 8) == [1, 0, 0]
    assert gray_bits(2, 8) == [1, 1, 0]  # gray_int(2) = 3
    assert gray_bits(3, 8) == [0, 1, 0]  # gray_int(3) = 2


def test_unary_bits():
    assert unary_bits(0, 4) == [1, 0, 0, 0]
    assert unary_bits(2, 4) == [0, 0, 1, 0]
    assert unary_bits(3, 4) == [0, 0, 0, 1]
    assert unary_bits(3, 7) == [0, 0, 0, 1, 0, 0, 0]

def test_unary_bits_invalid_level():
    with pytest.raises(ValueError):
        sb_bits(5, 4)
        sb_bits(5, 5)

def test_bitmask_subset_sb_and_gray():
    assert bitmask_subset(2, 8, "sb") == {0, 1, 2}
    assert bitmask_subset(2, 8, "gray") == {0, 1, 2}


def test_bitmask_subset_unary():
    assert bitmask_subset(0, 5, "unary") == {0}
    assert bitmask_subset(3, 5, "unary") == {3} 


def test_bitmask_subset_invalid_encoding():
    with pytest.raises(ValueError):
        bitmask_subset(1, 8, "banana")


def test_bits_for_level_sb():
    assert bits_for_level(5, 8, "sb") == [1, 0, 1]


def test_bits_for_level_gray():
    assert bits_for_level(2, 8, "gray") == [1, 1, 0]


def test_bits_for_level_unary():
    assert bits_for_level(3, 5, "unary") == [0, 0, 0, 1, 0]


def test_bits_for_level_invalid_encoding():
    with pytest.raises(ValueError):
        bits_for_level(1, 8, "banana")


def test_n_qubits():
    assert n_qubits(8, "sb") == 3
    assert n_qubits(8, "gray") == 3
    assert n_qubits(8, "unary") == 8


def test_n_qubits_invalid_encoding():
    with pytest.raises(ValueError):
        n_qubits(8, "banana")
