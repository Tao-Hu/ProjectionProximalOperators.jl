# Projection operator onto second order cone C = {(x,t) | |x|2 <= t}
# Usage: w = proj_soc(v)
# The first element of v is the t, and the rest elements are x
#
function proj_soc(v)

  # error check
  if v[1] < 0
    error("This is not a valid second order cone!");
  end

  # projection operator
  nx = norm(v[2:end]);
  n = length(v);
  if nx <= -v[1]
    w = zeros(n);
  elseif nx <= v[1]
    w = v;
  else
    coef = 0.5 * (1 + v[1] / nx);
    w = coef * [nx; v[2:end]];
  end

  return w;

end
