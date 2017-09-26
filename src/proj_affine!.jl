# in-place version of prox_affine
# Usage: prox_affine!(x, A, b)
#
function proj_affine!(x, A, b)

  BLAS.gemv!('N', 1.0, A, x, -1.0, b);
  BLAS.gemv!('N', -1.0, pinv(A), b, 1.0, x);
  return x;

end
