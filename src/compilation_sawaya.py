# Instead of qiskit compilation will try and implement according to description in the SI:
#
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

