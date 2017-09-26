# in-place version of prox_nuclear
# Usage: prox_nuclear!(X, lambda)
#
function prox_nuclear!(X, lambda)

  if lambda < 0

    return X;

  else

    (UX, evalX, VX) = svd(X);
    m = size(X, 1);
    n = size(X, 2);
    u = zeros(m);
    pu = pointer(u);
    v = zeros(n);
    pv = pointer(v);
    fill!(X, 0.0);

    for i = 1:min(m,n)
      tmpval = evalX - lambda;
      if tmpval > 0
        pU = pointer(UX) + (i - 1) * m * sizeof(Float64);
        BLAS.blascopy!(m, pU, 1, pu, 1);
        pV = pointer(VX) + (i - 1) * m * sizeof(Float64);
        BLAS.blascopy!(n, pV, 1, pv, 1);
        BLAS.ger!(tmpval, u, v, X);
      end
    end

    return X;

  end

end
