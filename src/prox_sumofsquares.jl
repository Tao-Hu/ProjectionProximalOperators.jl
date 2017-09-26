# Proximal operator of sum of squares g(x) = (lambda/2) * |x|_2^2
# Usage: y = prox_sumofsquares(x, lambda)
#
function prox_sumofsquares(x, lambda)

  y = x / (1 + lambda);
  return y;

end
