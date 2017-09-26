# Proximal operator of elastic net g(x) = lambda * (|x|_1 + (gamma/2)*|x|_2^2)
# Usage: y = prox_elasticnet(x, gamma, lambda)
#
function prox_elasticnet(x, gamma, lambda)

  y = similar(x);
  n = length(x);
  py = pointer(y);
  px = pointer(x);
  BLAS.blascopy!(n, px, 1, py, 1);
  coef = 1 / (1 + gamma * lambda);
  prox_l1!(y, lambda);
  BLAS.scal!(n, coef, y, 1);
  return y;

end
