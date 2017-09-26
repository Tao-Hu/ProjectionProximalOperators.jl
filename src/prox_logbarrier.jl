# Proximal operator of log barrier g(x) = -lambda * sum(log(x_i))
# Usage: y = prox_logbarrier(x, lambda)
#
function prox_logbarrier(x, lambda)

  y = similar(x);
  for i = 1 : length(x)
    y[i] = (x[i] + sqrt(x[i]^2 + 4 * lambda)) / 2;
  end

  return y;

end
