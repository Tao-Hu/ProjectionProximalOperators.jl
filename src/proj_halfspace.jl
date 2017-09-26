# Projection operator onto halfspace C = {x | a'x < b}
# Usage: y = proj_halfspace(x, a = avec, b = bval)
# Options:
#   a - vector a
#   b - scalar b
#
function proj_halfspace(x; a::Vector{Float64} = Float64[], b::Float64 = 0.0)

  # initialize
  n = length(x);
  if isempty(a)
    a = ones(n);
  end
  y = similar(x);
  px = pointer(x);
  py = pointer(y);
  BLAS.blascopy!(n, px, 1, py, 1);

  # pre calculate
  tmpval = BLAS.dot(n, a, 1, x, 1) - b;
  tmpval = max(tmpval, 0.0) / sumabs2(a);

  # projection operator
  BLAS.axpy!(-tmpval, a, y);

  return y;

end
