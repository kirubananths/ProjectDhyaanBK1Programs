function TE = TE_kraskov(X, Y, k_y, l_x, tau, kNN)

X = X(:); Y = Y(:);

% small noise (avoid ties)
X = X + 1e-10*randn(size(X));
Y = Y + 1e-10*randn(size(Y));

N = length(Y);

T = max(k_y, l_x)*tau + 1;
idx = T:N;
Npts = length(idx);

% future
Yf = Y(idx);

% past
Yp = zeros(Npts, k_y);
Xp = zeros(Npts, l_x);

for j = 1:k_y
    Yp(:,j) = Y(idx - j*tau);
end

for j = 1:l_x
    Xp(:,j) = X(idx - j*tau);
end

% spaces
Z_full = [Yf, Yp, Xp];
Z_yx   = [Yp, Xp];
Z_yy   = [Yf, Yp];
Z_y    = Yp;

% kNN search
[~, dist] = knnsearch(Z_full, Z_full, ...
    'K', kNN+1, 'Distance','chebychev');

eps = dist(:,end);

nx = zeros(Npts,1);
ny = zeros(Npts,1);
nz = zeros(Npts,1);

for i = 1:Npts
    
    nx(i) = sum(max(abs(Z_yy - Z_yy(i,:)),[],2) < eps(i)) - 1;
    ny(i) = sum(max(abs(Z_yx - Z_yx(i,:)),[],2) < eps(i)) - 1;
    nz(i) = sum(max(abs(Z_y  - Z_y(i,:)),[],2) < eps(i)) - 1;
    
end

TE = psi(kNN) + mean( psi(nz+1) - psi(ny+1) - psi(nx+1) );

end