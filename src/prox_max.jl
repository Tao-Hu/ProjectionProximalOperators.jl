# Proximal operator of max function g(x) = max(x)
# Usage: y = prox_max(x, lambda, tol = epsilon, nIters = n)
# Options:
#   tol - tolerence to stop iteration
#   nIters - maximum iterations
#
function prox_max(x, lambda; tol = 1e-8, nIters = 100)

  n = length(x);
  tl = minimum(x) - 1/n;
  tu = maximum(x);
  t0 = (tl + tu) / 2;

  for i = 1:nIters
    t0 = (tl + tu) / 2;
    tmpval1 = sum(max(x - t0, 0.0)) - lambda;
    tmpval2 = sum(max(x - tl, 0.0)) - lambda;
    if sign(tmpval1) == sign(tmpval2)
      tl = t0;
    else
      tu = t0;
    end
    if tu - tl <= tol
      break;
    end
  end

  return min(x, t0);

end
