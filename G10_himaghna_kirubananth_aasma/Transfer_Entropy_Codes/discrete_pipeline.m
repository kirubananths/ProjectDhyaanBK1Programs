function [TE, TE_corr, pval] = discrete_pipeline(X, Y, params, nSurr)

Xb = discretize_quantile(X, params.nBins);
Yb = discretize_quantile(Y, params.nBins);

[Yf, Yp, Xp] = build_states(Xb, Yb, params.k, params.l, params.tau);

Yp_idx = encode_state(Yp, params.nBins);
Xp_idx = encode_state(Xp, params.nBins);

TE = compute_TE_discrete(Yf, Yp_idx, Xp_idx);

TE_surr = zeros(nSurr,1);

for s = 1:nSurr
    shift = randi([round(0.1*length(Xp_idx)), round(0.9*length(Xp_idx))]);
    Xs = circshift(Xp_idx, shift);
    TE_surr(s) = compute_TE_discrete(Yf, Yp_idx, Xs);
end

TE_corr = TE - mean(TE_surr);
pval = mean(TE_surr >= TE);

end