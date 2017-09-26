# Algorithm is described in <http://icml2008.cs.helsinki.fi/papers/361.pdf>
# in-place version of proj_simplex
# Usage: proj_simplex!(x, b)
#
function proj_simplex!(x, b)

  # error check
  if b < 0
    error("Simplex radius can not be negative!");
  end

  # projection operator
  y = sort(x, rev = true);
  n = length(x);
  ysum = sum(y);
  rho = n;
  tmpval = 0.0;
  for i = n:-1:1
    tmpval = (b - ysum) / i;
    if y[i] + tmpval > 0
      rho = i;
      break;
    end
    ysum -= y[i];
  end
  lambda = tmpval;
  x = maximum(x + lambda, 0.0);

  return x;

end
