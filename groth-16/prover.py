# Groth-16 based prover for solution to x^3 + x - 30 = 0
# R1CS Constraints: 
# x * x = v1
# v1 * x = out - x + 30
#
# Transformation matrices satisfying Lw * Rw = Ow where w is the witness would be
# L = | 0 0 1 0 | R = | 0 0 1 0 | O = | 0   0  0  1 |  
#     | 0 0 0 1 |     | 0 0 1 0 |     | 30  1  -1 0 |
#
# The code given here very closely follows the notations used in the original Groth 16 paper

import numpy as np
from random import randint
import galois
from py_ecc.bn128 import add, multiply, curve_order, neg, pairing
from functools import reduce
from trusted_setup import get_secret_terms

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
u_polys = np.apply_along_axis(interpolate_column, 0, L_galois)
v_polys = np.apply_along_axis(interpolate_column, 0, R_galois)
w_polys = np.apply_along_axis(interpolate_column, 0, O_galois)


def generate_proof(sol):
  # Witness
  # Public: l = [0,1] and Private: l = [2, 3]

  a_public = GF(np.array([1, 0]))
  a_private = GF(np.array([sol, sol*sol]))
  a = GF(np.array([1, 0, sol, sol*sol]))

  def inner_product_polynomials_with_witness(polys, witness):
    mul_ = lambda x, y: x * y
    sum_ = lambda x, y: x + y
    return reduce(sum_, map(mul_, polys, witness))

  a_u = inner_product_polynomials_with_witness(u_polys, a)
  a_v = inner_product_polynomials_with_witness(v_polys, a)
  a_w = inner_product_polynomials_with_witness(w_polys, a)
  
  t = galois.Poly([1, curve_order - 1], field = GF) * galois.Poly([1, curve_order - 2], field = GF)
  h = (a_u * a_v - a_w) // t

  (sigma_1, sigma_2) = get_secret_terms(u_polys, v_polys, w_polys)
  (alpha_g1, beta_g1, delta_g1, x_g1_terms, private_internal_poly_eval_g1, public_internal_poly_eval_g1, powers_ht_g1) = sigma_1
  (beta_g2, gamma_g2, delta_g2, x_g2_terms) = sigma_2

  def poly_evaluation(vals, poly):
    mul_ = lambda x, y: multiply(x, y)
    sum_ = lambda x, y: add(x, y)
    return reduce(sum_, map(mul_, vals, poly))

  r = randint(0, curve_order)
  s = randint(0, curve_order)
  r_delta_g1 = multiply(delta_g1, r)
  s_delta_g1 = multiply(delta_g1, s)
  s_delta_g2 = multiply(delta_g2, s)

  A_g1 = add(alpha_g1, add(poly_evaluation(x_g1_terms, np.array(a_u.coeffs)), r_delta_g1))

  # Apparently a_v has degree 1 so a zero coeff must be inserted
  B_g2 = add(beta_g2, add(poly_evaluation(x_g2_terms, np.insert(np.array(a_v.coeffs), 0, 0)), s_delta_g2))

  # Required for C
  B_g1 = add(beta_g1, add(poly_evaluation(x_g1_terms, np.insert(np.array(a_v.coeffs), 0, 0)), s_delta_g1))

  h_t = poly_evaluation([powers_ht_g1], h.coeffs)

  # Terms to be added to create C
  A_g1_s = multiply(A_g1, s)
  B_g1_r = multiply(B_g1, r)
  rs_delta_g1 = multiply(delta_g1, r * s)
  C_sigma_term = add(poly_evaluation(private_internal_poly_eval_g1, a_private.base), h_t)

  C_g1 = add(C_sigma_term, add(A_g1_s, add(B_g1_r, neg(rs_delta_g1))))

  # Verification

  public_sigma_term = poly_evaluation(public_internal_poly_eval_g1, a_public.base)

  assert pairing(B_g2, A_g1) == pairing(beta_g2, alpha_g1) + pairing(gamma_g2, public_sigma_term) + pairing(delta_g2, C_g1)

  return (A_g1, B_g2, C_g1)

generate_proof(3)
  