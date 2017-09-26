# in-place version of prox_sumofsquares
# Usage: prox_sumofsquares!(x, lambda)
#
function prox_sumofsquares!(x, lambda)

  BLAS.scal!(length(x), 1/(1+lambda), x, 1);
  return x;

end
