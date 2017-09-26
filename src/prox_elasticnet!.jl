# in-place version of prox_elasticnet
# Usage: prox_elasticnet!(x, gamma, lambda)
#
function prox_elasticnet!(x, gamma, lambda)

  coef = 1 / (1 + gamma * lambda);
  prox_l1!(x, lambda);
  BLAS.scal!(n, coef, x, 1);
  return x;

end
