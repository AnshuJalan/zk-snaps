# QAP based verifier for X^3 + x - 30 = 0
# The prover sends encrypted polynomials AB = C
#
# The verifier checks for the pairing (A, B) to be equal to pairing(C, G2)

from py_ecc.bn128 import G2, pairing

# Accepts a proof in the form of (A, B, C)
def verify_proof(proof):
  (A, B, C) = proof
  return pairing(B, A) == pairing(G2, C)
