# QAP based prover for solution to x^3 + x - 30 = 0
# R1CS Constraints: 
# x * x = v1
# v1 * x = out - x + 30
#
# Transformation matrices satisfying Lw * Rw = Ow where w is the witness would be
# L = | 0 0 1 0 | R = | 0 0 1 0 | O = | 0   0  0  1 |  
#     | 0 0 0 1 |     | 0 0 1 0 |     | 30  1  -1 0 |
#

import numpy as np
import galois
from py_ecc.bn128 import add, multiply, curve_order
from functools import reduce
from trusted_setup import get_random_x_values

GF = galois.GF(curve_order)

L = np.array([[0, 0, 1, 0], [0, 0, 0, 1]])
R = np.array([[0, 0, 1, 0], [0, 0, 1, 0]])
O = np.array([[0, 0, 0, 1], [30, 1, curve_order - 1, 0]])

# Field matrices
L_galois = GF(L)
R_galois = GF(R)
O_galois = GF(O)

def interpolate_column(col):
  xs = GF(np.array([1,2]))
  return galois.lagrange_poly(xs, col)

# L, R, O resolved to polynomials
U_polys = np.apply_along_axis(interpolate_column, 0, L_galois)
V_polys = np.apply_along_axis(interpolate_column, 0, R_galois)
W_polys = np.apply_along_axis(interpolate_column, 0, O_galois)


def generate_proof(x):
  # solution vector
  witness = GF(np.array([1, 0, x, x * x]))

  def inner_product_polynomials_with_witness(polys, witness):
    mul_ = lambda x, y: x * y
    sum_ = lambda x, y: x + y
    return reduce(sum_, map(mul_, polys, witness))

  term_1 = inner_product_polynomials_with_witness(U_polys, witness)
  term_2 = inner_product_polynomials_with_witness(V_polys, witness)
  term_3 = inner_product_polynomials_with_witness(W_polys, witness)

  t = galois.Poly([1, curve_order - 1], field = GF) * galois.Poly([1, curve_order - 2], field = GF)
  h = (term_1 * term_2 - term_3) // t

  evaluation_target, evalutation_vals = get_random_x_values()

  def poly_evaluation(vals, poly):
    mul_ = lambda x, y: multiply(x, y)
    sum_ = lambda x, y: add(x, y)
    return reduce(sum_, map(mul_, vals, poly))

  A = poly_evaluation(evalutation_vals[0], np.array(term_1.coeffs))

  # Apparently term_2 has degree 1 so a zero coeff must be inserted
  B = poly_evaluation(evalutation_vals[1], np.insert(np.array(term_2.coeffs), 0, 0))

  t_eval = t(evaluation_target)
  h_t = poly_evaluation([multiply(evalutation_vals[0][i], np.array(t_eval.base)[0]) for i in range(2)], h.coeffs)

  C = add(poly_evaluation(evalutation_vals[0], np.array(term_3.coeffs)), h_t)

  return (A, B, C)
  