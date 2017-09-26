# in-place version of prox_max
# Usage: prox_max!(x, lambda, tol = epsilon, nIters = n)
#
function prox_max!(x, lambda; tol = 1e-8, nIters = 100)

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

  for i = 1:n
    x[i] = min(x[i], t0);
  end
  return x;

end
