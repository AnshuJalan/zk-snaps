pragma circom 2.1.6;

// Circuit for the polynomial x^3 + x - 30 = 0
template Polynomial() {
  signal input x;
  
  signal v1; 
  signal output out;

  v1 <== x * x;
  out <== (v1 * x) + x - 30;
}

component main = Polynomial();