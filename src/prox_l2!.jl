# in-place versionof prox_l2
# Usage: prox_l2!(x, lambda, grpSize)
#
function prox_l2!(x, lambda, grpSize)

  g = length(grpSize);
  px = pointer(x);

  for i = 1:g
    coef = max(1 - lambda[i] / BLAS.nrm2(grpSize[i], px, 1), 0.0);
    BLAS.scal!(grpSize[i], coef, px, 1);
    px += grpSize[i] * sizeof(Float64);
  end

  return x;

end
