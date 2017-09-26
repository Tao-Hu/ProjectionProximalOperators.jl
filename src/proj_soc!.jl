# in-place version of proj_soc
# Usage: proj_soc!(v)
#
function proj_soc!(v)

  # error check
  if v[1] < 0
    error("This is not a valid second order cone!");
  end

  # projection operator
  nx = norm(v[2:end]);
  n = length(v);
  if nx <= -v[1]
    fill!(v, 0.0);
  elseif nx >= abs(v[1])
    coef = 0.5 * (1 + v[1] / nx);
    v[1] = nx;
    BLAS.scal!(n, coef, v, 1);
  end

  return v;

end
