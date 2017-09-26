# Proximal operator of l1 norm g(x) = lambda * |x|_1
# Usage: y = prox_l1(x, lambda)
#
function prox_l1(x, lambda)

  y = similar(x);

  for i = 1:length(x)
    y[i] = copysign(max(abs(x[i]) - lambda, 0.0), x[i]);
  end

  return y;

end
