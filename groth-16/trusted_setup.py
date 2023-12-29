# Trusted setup for Groth 16
#
# Generates σ_{1} and σ_{2} as defined in section 3.2 of groth16 paper

import galois 
import numpy as np
from random import randint
from py_ecc.bn128 import multiply, curve_order, G1, G2

GF = galois.GF(curve_order)

# Public target 
t = galois.Poly([1, curve_order - 1], field = GF) * galois.Poly([1, curve_order - 2], field = GF)

# Accepts polynomials u, v, w same as defined in the paper
def get_secret_terms(u_polys, v_polys, w_polys):
  # Safer to use a stronger source of randomness or pre-select hidden random numbers
  alpha = GF(randint(0, curve_order))
  beta = GF(randint(0, curve_order))
  delta = GF(randint(0, curve_order))
  gamma = GF(randint(0, curve_order))
  x = GF(randint(0, curve_order))

  def eval_and_sum(index, denom):
    return ((beta * u_polys[index](x)) + (alpha * v_polys[index](x)) + w_polys[index](x)) / denom
  
  # Evaluation of polys within sigma of private and public input
  private_internal_poly_eval = [eval_and_sum(i, gamma) for i in [0, 1]]
  public_internal_poly_eval = [eval_and_sum(i, delta) for i in [2, 3]]

  powers_ht = t(x) / delta

  # Takes a field element and converts it to raw integer
  def extract(n):
    return np.array(n.flatten())[0]

  alpha_g1 = multiply(G1, extract(alpha))
  beta_g1 = multiply(G1, extract(beta))
  delta_g1 = multiply(G1, extract(delta))
  x_g1 = multiply(G1, extract(x))
  private_internal_poly_eval_g1 = [multiply(G1, extract(n)) for n in private_internal_poly_eval]
  public_internal_poly_eval_g1 = [multiply(G1, extract(n)) for n in public_internal_poly_eval]
  power_ht_g1 = multiply(G1, extract(powers_ht))

  # This can of course be made efficient by caching extracted values, but that is not the 
  # goal here for the demonstration
  beta_g2 = multiply(G2, extract(beta))
  gamma_g2 = multiply(G2, extract(gamma))
  delta_g2 = multiply(G2, extract(delta))
  x_g2 = multiply(G2, extract(x))

  sigma_1 = (alpha_g1, beta_g1, delta_g1, [x_g1, G1], private_internal_poly_eval_g1, public_internal_poly_eval_g1, power_ht_g1)
  sigma_2 = (beta_g2, gamma_g2, delta_g2, [x_g2, G2])

  return (sigma_1, sigma_2)