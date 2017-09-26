# in-place version of prox_l1
# Usage: prox_l1!(x, lambda)
#
function prox_l1!(x, lambda)

  for i = 1:length(x)
    x[i] = copysign(max(abs(x[i]) - lambda, 0.0), x[i]);
  end

  return x;

end
