%% =========================================
%   Transfer Entropy Comparison Script
% =========================================

%clear; clc;
disp('Started');
%% -------------------------
% 1. Simulate signals
% -------------------------
N = 2000;

Z = randn(N,1);  % common source
X = filter(1, [1 -0.7], Z) + 0.3*randn(N,1);
Y = filter(1, [1 -0.7], Z) + 0.3*randn(N,1);

% Add real causal influence
for t = 2:N
    Y(t) = Y(t) + 0.4*X(t-1);
end
disp('Started');
%% -------------------------
% 2. Parameters
% -------------------------
params.nBins = 4;
params.k = 2;
params.l = 2;
params.tau = 1;

kNN = 4;
nSurr = 50;
disp('Started');
%% -------------------------
% 3. DISCRETE TE
% -------------------------
[TE_d, TEc_d, p_d] = discrete_pipeline(X, X, params, nSurr);
disp('Done dis');
%% -------------------------
% 4. KRASKOV TE
% -------------------------
[TE_k, TEc_k, p_k] = kraskov_pipeline(X, Y, params, kNN, nSurr);
disp('Done kr');
%% -------------------------
% 5. CONDITIONAL KRASKOV TE
% -------------------------
[TE_ck, TEc_ck, p_ck] = conditional_pipeline(X, Y, Z, params, kNN, nSurr);
disp('Done ckr');
%% -------------------------
% 6. DISPLAY
% -------------------------
fprintf('\n===== RESULTS =====\n');

fprintf('\n--- Discrete TE ---\n');
fprintf('TE      = %.4f\n', TE_d);
fprintf('TE_corr = %.4f\n', TEc_d);
fprintf('p       = %.4f\n', p_d);

fprintf('\n--- Kraskov TE ---\n');
fprintf('TE      = %.4f\n', TE_k);
fprintf('TE_corr = %.4f\n', TEc_k);
fprintf('p       = %.4f\n', p_k);

fprintf('\n--- Conditional Kraskov TE ---\n');
fprintf('TE      = %.4f\n', TE_ck);
fprintf('TE_corr = %.4f\n', TEc_ck);
fprintf('p       = %.4f\n', p_ck);