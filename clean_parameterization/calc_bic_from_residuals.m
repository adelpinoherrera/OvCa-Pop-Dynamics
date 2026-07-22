function BIC = calc_bic_from_residuals(resvec, k)
    N = numel(resvec);
    SSE = sum(resvec.^2);
    BIC = N * log(SSE / N) + k * log(N);
end