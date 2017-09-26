# in-place version of prox_halfspace
# Usage: proj_halfspace!(x, a = avec, b = bval)
#
function proj_halfspace!(x; a::Vector{Float64} = Float64[], b::Float64 = 0.0)

  # initialize
  n = length(x);
  if isempty(a)
    a = ones(n);
  end

  # pre calculate
  tmpval = BLAS.dot(n, a, 1, x, 1) - b;
  tmpval = max(tmpval, 0.0) / sumabs2(a);

  # projection operator
  BLAS.axpy!(-tmpval, a, x);

  return x;

end
