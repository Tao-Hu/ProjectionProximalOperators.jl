# Projection operator onto positive semidefinite cone
# Usage: Y = proj_psdc(X)
#
function proj_psdc(X)

  (d, U) = eig(X);
  n = size(X, 1);
  Y = zeros(n, n);
  u = zeros(n);
  pu = pointer(u);

  for i = 1:n
    if d[i] > 0
      pU = pointer(U) + (i - 1) * n * sizeof(Float64);
      BLAS.blascopy!(n, pU, 1, pu, 1);
      BLAS.ger!(d[i], u, u, Y);
    end
  end

  return Y;

end
