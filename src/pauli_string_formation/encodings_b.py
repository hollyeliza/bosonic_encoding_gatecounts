import math

def n_qubits_sb(d: int) -> int:
    """"
    Returns the number of qubits needed to encode Fock states up to d
    (d = N_ max + 1) in standard binary.
    
    It is N_ max + 1 because we need |0> to |N_max>.
    
    If we want to encode up to |7> (|N_ max>) we need 8 bitstrings: 
    
    |0>, |1>, |2>, |3>, |4>, |5>, |6>, |7>, |8>
    -->
    |000>, |001>, |010>, |011>, |100>, |101>, |110>,|111>
    """
    return math.ceil(math.log2(d))


def sb_bits(l: int, d: int) -> list[int]:
    """
    Returns the represention of each integer l in standard binary,
    for a given number of qubits. This returns the bits of l as a list, read 
    from right to left (least significant bit (0) first).
    MSB (most signicant bit) typically first but here LSB is first.

    Raises
    ------
    ValueError
        If d <= 0 or if l < d
    """
    if l >= d:
        raise ValueError(f"l must satisfy l < d, but got l={l}, d={d}")
    
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
    Returns the representation of each integer l in its Gray code representation 
    (as an integer). In Gray code, adjacent values differ by only one bit.
    """
    return l ^ (l >> 1) # The ^ is XOR

# print(gray_int(8)) # 8 in gray is 12 in binary (so machines stores gray encoding of 8 as 12)

def gray_bits(l: int, d: int) -> list[int]:
    """"
    Returns the represention of each integer l in the gray encoding as a list with MSB last.
    """
    g = gray_int(l)
    nq = n_qubits_sb(d)
    return [(g >> i) & 1 for i in range(nq)]


def unary_bits(l: int, d: int) -> list[int]:
    """"
    Returns the represention of each integer l in the unary encoding as a list with MSB last
    (read from right to left). This is easy and just makes a list of 0's and shoves a 1 at position l.
    This returns the bits of l in the unary encoding as a list, read from right to left.

    Note that for unary encoding a lot more bits are needed (d = N_max + 1).
    """

    if l >= d:
        raise ValueError(f"l must satisfy l < d, but got l={l}, d={d}")
    
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
        return {l} # Remember thought that you read these strings right to left. 
    else:
        raise ValueError(f"Unknown encoding: {encoding}")
    
# print(f'bitstring sb: {sb_bits(5, 4)}; bitmask subset for sb example: {bitmask_subset(5,4, "sb")}')
# print(f'bitstring gray {gray_bits(5, 4)}; bitmask subset for gray example: {bitmask_subset(5,4, "gray")}')
# print(f'bitstring unary: {unary_bits(5, 6)}; bitmask subset for unary example: {bitmask_subset(5,6, "unary")}')

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
    
