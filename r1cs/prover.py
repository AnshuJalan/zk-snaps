# R1CS based prover for solution to x^3 + x - 30 = 0
# Constraints: 
# v1 = x * x
# out - x + 30 = v1 * x

# The prover sends EC encrypted points in the witness to the verifier
# witness: [1, out, x, v1] where out = 0

from py_ecc.bn128 import G1, G2, multiply

# Returns an array of variables of two witnesses
# At index 0 -> wG1 (without 1 and out)
# At index 1 -> wG2 ( " )
def generate_proof(x):
    # Intermediary
    v1 = x * x

    # Encrypt the variables
    v1_G1 = multiply(G1, v1)
    v1_G2 = multiply(G2, v1)
    x_G1 = multiply(G1, x)
    x_G2 = multiply(G2, x)

    return [[x_G1, v1_G1], [x_G2, v1_G2]]
