# Crafts a random number between 0 and field_modulus

from py_ecc.bn128 import G1, G2, multiply, field_modulus
from random import randint

# Based on the R1CS of x^3 + x - 30 which requires degree 1 polynomial to be represented,
# this function returns <G1, r.G1>, <G2, r.G2>
def get_random_x_values():
  r = 5
  return r, [[multiply(G1, r), G1], [multiply(G2, r), G2]]