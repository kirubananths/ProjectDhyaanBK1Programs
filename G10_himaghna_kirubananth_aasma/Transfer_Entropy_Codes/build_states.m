function [Yf, Yp, Xp] = build_states(Xb, Yb, k, l, tau)

T = max(k,l)*tau + 1;
idx = T:length(Yb);

Yf = Yb(idx);

Yp = zeros(length(idx), k);
Xp = zeros(length(idx), l);

for j = 1:k
    Yp(:,j) = Yb(idx - j*tau);
end

for j = 1:l
    Xp(:,j) = Xb(idx - j*tau);
end

end