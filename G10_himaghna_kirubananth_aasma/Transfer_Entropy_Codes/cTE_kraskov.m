function TE = cTE_kraskov(X, Y, Z, k_y, l_x, l_z, tau, kNN)

% small noise (avoid ties)
X = X(:) + 1e-10*randn(size(X(:)));
Y = Y(:) + 1e-10*randn(size(Y(:)));
Z = Z(:) + 1e-10*randn(size(Z(:)));

N = length(Y);

T = max([k_y, l_x, l_z])*tau + 1;
idx = T:N;

Npts = length(idx);

% --- Build embeddings ---
Yf = Y(idx);

Yp = zeros(Npts, k_y);
Xp = zeros(Npts, l_x);
Zp = zeros(Npts, l_z);

for j = 1:k_y
    Yp(:,j) = Y(idx - j*tau);
end

for j = 1:l_x
    Xp(:,j) = X(idx - j*tau);
end

for j = 1:l_z
    Zp(:,j) = Z(idx - j*tau);
end

% --- Define spaces ---
Z_full = [Yf, Yp, Xp, Zp];   % full
A = [Yp, Zp];                % conditioning
B = [Yp, Xp, Zp];            % without Yf
C = [Yf, Yp, Zp];            % without Xp

% --- kNN in full space ---
[idx_knn, dist] = knnsearch(Z_full, Z_full, ...
    'K', kNN+1, 'Distance','chebychev');

eps = dist(:,end);

% --- Count neighbors ---
nA = zeros(Npts,1);
nB = zeros(Npts,1);
nC = zeros(Npts,1);

for i = 1:Npts
    
    nA(i) = sum(max(abs(A - A(i,:)),[],2) < eps(i)) - 1;
    nB(i) = sum(max(abs(B - B(i,:)),[],2) < eps(i)) - 1;
    nC(i) = sum(max(abs(C - C(i,:)),[],2) < eps(i)) - 1;
    
end

% --- TE estimate ---
TE = psi(kNN) + mean( psi(nA+1) - psi(nB+1) - psi(nC+1) );

end