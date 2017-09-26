# in-place version of prox_affine
# Usage: prox_affine!(x, b, lambda)
#
function prox_affine!(x, b, lambda)

  x = x - lambda * b;
  return x;

end
