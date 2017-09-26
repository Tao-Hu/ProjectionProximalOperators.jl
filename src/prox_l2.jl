# Proximal operator of l2 norm g(x) = lambda * (sum of |x_g|_2)
# Usage: y = prox_l2(x, lambda, grpSize)
# Options:
#   grpSize - size for each group
#
function prox_l2(x, lambda, grpSize)

  g = length(grpSize);
  y = similar(x);
  py = pointer(y);
  px = pointer(x);
  BLAS.blascopy!(length(x), px, 1, py, 1);

  for i = 1:g
    coef = max(1 - lambda[i] / BLAS.nrm2(grpSize[i], px, 1), 0.0);
    BLAS.scal!(grpSize[i], coef, py, 1);
    px += grpSize[i] * sizeof(Float64);
    py += grpSize[i] * sizeof(Float64);
  end

  return y;

end
