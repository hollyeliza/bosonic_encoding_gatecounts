import math

def n_qubits_sb(d: int) -> int:
    """"
    For a cut-off of d, how many qubits are needed to encode these Fock 
    states in standard binary.
    """
    return math.ceil(math.log2(d))

def sb_bits(l: int, d: int) -> list[int]:
    """
    What is the representation of each integer l in standard binary,
    for a given number of qubits.
    This returns the bits of l as a list, read from right to left
    (least significant bit (0) first).
    """
    nq = n_qubits_sb(d)
    bitstring = []
    for i in range(nq):
        # Extract the i-th bit of l (little-endian order)
        last_bit = (l >> i) & 1 # takes the last bit of l and comapares it to 1 (l stored in binry)
        bitstring.append(last_bit) # if 1 then 1 is appended; if 0 then 0 appended
        # print(f'i = {i}, extracted bit = {last_bit}, bitstring so far = {bitstring}')
    return bitstring # you read bits from right to left!

def gray_int(l: int) -> int:
    """
    Convert integer l into its Gray code representation (as an integer).
    In Gray code, adjacent values differ by only one bit.
    """
    return l ^ (l >> 1) # The ^ is XOR

# print(gray_int(8)) # 8 in gray is 12 in binary (so machines stores gray encoding of 8 as 12)

def gray_bits(l: int, d: int) -> list[int]:
    """"
    This returns the bits of l in the gray encoding as a list, read from right to left.
    """
    g = gray_int(l)
    nq = n_qubits_sb(d)
    return [(g >> i) & 1 for i in range(nq)]

def unary_bits(l: int, d: int) -> list[int]:
    """"
    Easy - makes a list of 0's and shoves a 1 at position l.
    This returns the bits of l in the unary encoding as a list, read from right to left.
    Note that for unary encoding a lot more bits are needed - n_qubits_sb not used. d is l + 1
    """
    bits = [0] * d
    bits[l] = 1
    return bits

def bitmask_subset(l: int, d: int, encoding: str) -> set[int]:
    """
    Returns the labelling of the qubits associated with the encoding.
    For more compact encodings (sba nd gray) it is all of them - 
    all the qubits are needed to know what level this encoding represents.
    For the unary encoding only the qubit corresponding to the only 1 is important.
    """
    if encoding in ("sb", "gray"):
        nq = n_qubits_sb(d)
        return set(range(nq))
    elif encoding == "unary":
        return {l}
    else:
        raise ValueError(f"Unknown encoding: {encoding}")


def bits_for_level(l: int, d: int, encoding: str) -> list[int]:
    """"
    Returns the bitstring for the l (occupation number) in the correct encoding.
    Read from right to left!!
    Uses all of the above helper functions.
    """
    if encoding == "sb":
        return sb_bits(l, d)
    elif encoding == "gray":
        return gray_bits(l, d)
    elif encoding == "unary": 
        return unary_bits(l, d)
    else:
        raise ValueError(f"Unknown encoding: {encoding}")

# print(bits_for_level(5, 6, "unary")) # In reality the 1 is at the lhs (read from right to left!)

def n_qubits(d: int, encoding: str) -> int:
    """
    Returns the number of qubits needed to represent each of the encodings.
    """
    if encoding in ("sb", "gray"):
        return n_qubits_sb(d)
    elif encoding == "unary":
        return d
    else:
        raise ValueError(f"Unknown encoding: {encoding}")
    
