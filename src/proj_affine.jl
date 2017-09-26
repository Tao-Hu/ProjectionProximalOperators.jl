# Projection operator onto affine set C = {x | Ax = b}
# Usage: y = proj_affine(x, A, b)
#
function proj_affine(x, A, b)

  (m, n) = size(A);
  tmpvec = similar(b);
  pb = pointer(b);
  pvec = pointer(tmpvec);
  BLAS.blascopy!(m, pb, 1, pvec, 1);
  BLAS.gemv!('N', 1.0, A, x, -1.0, tmpvec);
  y = x - BLAS.gemv('N', pinv(A), tmpvec);
  return y;

end
