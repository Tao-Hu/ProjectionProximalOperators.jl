# Proximal operator of quadratic function g(x) = lambda * (0.5*x'*A*x + b'*x + c)
# Usage: y = prox_quad(x, A, b, alg = algorithm, tolX = epsilon, xInit = x0, nIter = n)
# Options:
#   alg - algorithm to solve this problem. For small scale problem, you can specify
#         it as "QR"; for large scale problem, "CG" is prefered. Default value is "QR"
#   tolX - tolerence to stop iteration when alg is specified as "CG"
#   xInit - initial value for iteration when alg is specified as "CG"
#   nIter - maximum iterations
#
function prox_quad(x, A, b, lambda; alg::String = "QR", tolX::Float64 = 1e-5,
                   xInit::Vector{Float64} = Float64[], nIter::Int = 0)

  # error check
  n = length(x);
  if size(A, 1) != n || size(A, 2) != n
    error("Dimension of matrix A must match with x!");
  end
  if length(b) != n
    error("Dimension of vector b must match with x!");
  end

  # proximal operator
  tmpA = eye(n) + lambda * A;
  tmpb = x - lambda * b;
  if alg == "QR"

    y = tmpA \ tmpb;
    return y;

  elseif alg == "CG"

    if isempty(xInit)
      xInit = ones(n);
    end
    if nIter == 0
      nIter = n;
    end
    y = copy(xInit);
    r = similar(b);
    pr = pointer(r);
    pb = pointer(b);
    BLAS.blascopy!(n, pb, 1, pr, 1);
    BLAS.gemv!('N', 1.0, tmpA, xInit, -1.0, r);
    rsold = sumabs2(r);
    p = BLAS.scal(n, -1.0, r, 1);
    Ap = BLAS.gemv('N', tmpA, p);

    for i = 1:nIter
      alpha = rsold / BLAS.dot(n, p, 1, Ap, 1);
      y += alpha * p;
      r += alpha *Ap;
      rsnew = sumabs2(r);
      if sqrt(rsnew) <= tolX
        break;
      end
      beta = rsnew / rsold;
      p = beta * p - r;
      rsold = rsnew;
      BLAS.gemv!('N', 1.0, tmpA, p, 0.0, Ap);
    end

    return y;

  else

    error("Unsupported algorithm to solve linear equation system.");

  end

end
