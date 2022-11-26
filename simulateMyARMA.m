function y = simulateMyARMA(N, sigma2, ARMA_poly)
    rng(0);
    e = sqrt(sigma2) * randn(N,1);
    y = filter(ARMA_poly.c, ARMA_poly.a, e);
end