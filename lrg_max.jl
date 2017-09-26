function lrgmaxpath(y, X, λList; alg::String = "FISTA", nIters::Int = 1000,
                 tol::Float64 = 1e-4)

  p = size(X, 2);
  nλ = length(λList);
  βhat = zeros(p, nλ);

  # pre-compute
  XtX = BLAS.gemm('T', 'N', X, X);
  Xty = BLAS.gemv('T', X, y);
  stepSize = 1 / eigmax(XtX);

  for i = nλ:-1:1
    λ = λList[i];
    if i == nλ
      βhat[:, i] = lrgmax(y, X, λ, alg = alg, nIters = nIters, tol = tol,
                          stepSize = stepSize, Xty = Xty, XtX = XtX);
    else
      βhat[:, i] = lrgmax(y, X, λ, alg = alg, nIters = nIters, tol = tol,
                          stepSize = stepSize, Xty = Xty, XtX = XtX,
                          β0 = βhat[:, i+1]);
    end
  end

  return βhat;

end

function lrgmax(y, X, λ; alg::String = "FISTA", nIters::Int = 1000,
                tol::Float64 = 1e-4, β0::Vector{Float64} = zeros(size(X, 2)),
                stepSize::Float64 = 1.0, Xty::Vector{Float64} = Float64[],
                XtX::Matrix{Float64} = [Float64[] Float64[]])

  (n, p) = size(X);

  if alg == "Convex.jl"

    βhat = Convex.Variable(size(X, 2));
    problem = Convex.minimize(0.5 * Convex.sum_squares(y - X * βhat)
                              + λ * Convex.maximum(βhat));
    Convex.solve!(problem);
    return βhat.value;

  elseif alg == "FISTA"

    # Initial
    s = stepSize;
    βhat = copy(β0);
    pβhat = pointer(βhat);
    residual = copy(y);
    BLAS.gemv!('N', -1.0, X, βhat, 1.0, residual);
    py = pointer(y);
    presidual = pointer(residual);
    obj = 0.5 * sumabs2(residual) + λ * maximum(βhat);
    p1b = copy(βhat);
    p2b = copy(βhat);
    pp1b = pointer(p1b);
    pp2b = pointer(p2b);
    pXty = pointer(Xty);

    # FISTA algorithm
    for i = 1:nIters

      # extrapolation
      oldobj = obj;
      coef1 = (2 * i - 1) / (i + 1);
      coef2 = (2 - i) / (i + 1);
      tmpy = coef1 * p1b + coef2 * p2b;
      BLAS.blascopy!(p, pXty, 1, pβhat, 1);
      BLAS.gemv!('N', -1.0, XtX, tmpy, 1.0, βhat);
      BLAS.scal!(p, s, βhat, 1);
      BLAS.axpy!(p, 1.0, tmpy, 1, βhat, 1);

      # update βhat
      prox_max!(βhat, s*λ);

      # check convergence
      BLAS.blascopy!(n, py, 1, presidual, 1);
      BLAS.gemv!('N', -1.0, X, βhat, 1.0, residual);
      obj = 0.5 * sumabs2(residual) + λ * maximum(βhat);
      if abs(oldobj - obj) <= tol * (oldobj + 1)
        break
      else
        BLAS.blascopy!(p, pp1b, 1, pp2b, 1);
        BLAS.blascopy!(p, pβhat, 1, pp1b, 1);
      end
    end

    return βhat;

  else

    error("Unsupported algorithm!");

  end

end