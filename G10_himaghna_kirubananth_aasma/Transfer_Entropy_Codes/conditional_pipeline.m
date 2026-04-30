function [TE, TE_corr, pval] = conditional_pipeline(X, Y, Z, params, kNN, nSurr)

TE = cTE_kraskov(X, Y, Z, params.k, params.l, params.l, params.tau, kNN);

TE_surr = zeros(nSurr,1);

for s = 1:nSurr
    shift = randi([round(0.1*length(X)), round(0.9*length(X))]);
    Xs = circshift(X, shift);
    TE_surr(s) = cTE_kraskov(Xs, Y, Z, params.k, params.l, params.l, params.tau, kNN);
end

TE_corr = TE - mean(TE_surr);
pval = mean(TE_surr >= TE);

end