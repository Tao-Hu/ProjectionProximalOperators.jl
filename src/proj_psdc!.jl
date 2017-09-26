# in-place version of proj_psdc
# Usage: proj_psdc!(X)
#
function proj_psdc!(X)

  (d, U) = eig(X);
  n = size(X, 1);
  fill!(X, 0.0);
  u = zeros(n);
  pu = pointer(u);

  for i = 1:n
    if d[i] > 0
      pU = pointer(U) + (i - 1) * n * sizeof(Float64);
      BLAS.blascopy!(n, pU, 1, pu, 1);
      BLAS.ger!(d[i], u, u, X);
    end
  end

  return X;

end
