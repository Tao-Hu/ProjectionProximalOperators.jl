# in-place version of prox_logbarrier
# Usage: prox_logbarrier!(x, lambda)
#
function prox_logbarrier!(x, lambda)

  for i = 1 : length(x)
    x[i] = (x[i] + sqrt(x[i]^2 + 4 * lambda)) / 2;
  end

  return x;

end
