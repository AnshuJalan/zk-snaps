# R1CS based verifier for X^3 + x - 30 = 0
# The prover sends a partly encrypted witness [1, out, x, v1]
# 
# If,
#           v1 = x * x
# out - x + 30 = v1 * x
#
# Transformation matrices satisfying Aw * Bw = Cw would be
# A = | 0 0 1 0 | B = | 0 0 1 0 | C = | 0   0  0  1 |  
#     | 0 0 0 1 |     | 0 0 1 0 |     | 30  1  -1 0 |
#

from py_ecc.bn128 import G1, G2, pairing, add, multiply, neg

A = [[0, 0, 1, 0], [0, 0, 0, 1]]
B = [[0, 0, 1, 0], [0, 0, 1, 0]]
C = [[0, 0, 0, 1], [30, 1, -1, 0]]

out_G1 = multiply(G1, 0)
out_G2 = multiply(G2, 0)
one_G1 = G1
one_G2 = G2

# Tailor made dot product for the matrices in the verifier
def dot_product(m1, m2):
    res = []
    for i in range(2):
        s = multiply(m2[0], m1[i][0])
        for j in range(1, 4):
            witness_element = neg(m2[j]) if m1[i][j] < 0 else m2[j]
            s = add(s, multiply(witness_element, abs(m1[i][j])))
        res.append(s)
    return res

# Accepts a proof in the form of an R1CS witness [[x_G1, v1_G1], [x_G2, v1_G2]]
def verify_proof(proof):
    assert len(proof) == 2
    assert (len(proof[0]) == 2) and (len(proof[1]) == 2)

    # Fully encrypted proofs (witness)
    w_G1 = [one_G1, out_G1, proof[0][0], proof[0][1]]
    w_G2 = [one_G2, out_G2, proof[1][0], proof[1][1]]

    # Must be of the same orders (1 x 2)
    Aw_G1 = dot_product(A, w_G1)
    Bw_G2 = dot_product(B, w_G2)
    Cw_G1 = dot_product(C, w_G1)

    # Pairing applied to individual elements of Aw_G1 and Bw_G2
    # Based on the constraints above, the length must be 2
    lhs = [pairing(Bw_G2[i], Aw_G1[i]) for i in range(2)]
    rhs = [pairing(G2, Cw_G1[i]) for i in range(2)]

    return lhs == rhs
