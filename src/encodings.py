import math


def n_qubits_sb(d: int) -> int:
    return math.ceil(math.log2(d))


def sb_bits(l: int, d: int) -> list[int]:
    nq = n_qubits_sb(d)
    return [(l >> i) & 1 for i in range(nq)]


def gray_int(l: int) -> int:
    return l ^ (l >> 1)


def gray_bits(l: int, d: int) -> list[int]:
    g = gray_int(l)
    nq = n_qubits_sb(d)
    return [(g >> i) & 1 for i in range(nq)]


def unary_bits(l: int, d: int) -> list[int]:
    bits = [0] * d
    bits[l] = 1
    return bits


def bitmask_subset(l: int, d: int, encoding: str) -> set[int]:
    if encoding in ("sb", "gray"):
        nq = n_qubits_sb(d)
        return set(range(nq))
    elif encoding == "unary":
        return {l}
    else:
        raise ValueError(f"Unknown encoding: {encoding}")


def bits_for_level(l: int, d: int, encoding: str) -> list[int]:
    if encoding == "sb":
        return sb_bits(l, d)
    elif encoding == "gray":
        return gray_bits(l, d)
    elif encoding == "unary":
        return unary_bits(l, d)
    else:
        raise ValueError(f"Unknown encoding: {encoding}")


def n_qubits(d: int, encoding: str) -> int:
    if encoding in ("sb", "gray"):
        return n_qubits_sb(d)
    elif encoding == "unary":
        return d
    else:
        raise ValueError(f"Unknown encoding: {encoding}")
    

gray = n_qubits(32, "gray")
print(f'The number of qubits needed for a bosonic cut off 32 is {gray}')