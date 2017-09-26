# Projection operator onto box set C = {x | l <= x <= u}
# Usage: y = proj_box(x, lb = lbvec, ub = ubvec)
# Options:
#   lb - vector of lower bound
#   ub - vector of upper bound
#
function proj_box(x; lb::Vector{Float64} = Float64[],
                  ub::Vector{Float64} = Float64[])

  n = length(x);
  y = similar(x);

  # initial lower bound and upper bound
  if isempty(lb)
    lb = zeros(n);
    fill!(lb, -Inf);
  end
  if isempty(ub)
    ub = zeros(n);
    fill!(ub, Inf);
  end

  # error check
  if length(lb) != n
    error("proj.box: dimension of lower bound vector does not match!");
  end
  if length(ub) != n
    error("proj.box: dimension of upper bound vector does not match!");
  end

  # projection operator
  for i = 1 : n
    if lb[i] > ub[i]
      error("proj.box: lower bound can not be larger than upper bound!");
    end
    y[i] = max(lb[i], min(x[i], ub[i]));
  end

  return y;
end
