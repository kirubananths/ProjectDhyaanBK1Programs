function Xb = discretize_quantile(X, nBins)
%edges = quantile(X, linspace(0,1,nBins+1));
%edges(1) = -inf; edges(end) = inf;
Xb = discretize(X, nBins);
end