# Proximal operator of affine function g(x) = lambda * (b'x + c)
# Usage: y = prox_affine(x, b, lambda)
#
function prox_affine(x, b, lambda)

  y = x - lambda * b;
  return y;

end
