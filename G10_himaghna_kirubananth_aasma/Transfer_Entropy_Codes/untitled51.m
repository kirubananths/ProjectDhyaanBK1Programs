alpha = 0.05;

%% =========================
% 1. THRESHOLD EACH TRIAL
%% =========================

[nElec, ~, nTrials] = size(mCo1b);

TE_sig = zeros(size(mCo1b));

for t = 1:nTrials
    mask = mpv1b(:,:,t) < alpha;
    TE_sig(:,:,t) = mCo1b(:,:,t) .* mask;
end

%% =========================
% 2. AVERAGE ONLY SIGNIFICANT VALUES
%% =========================

TE_sum = sum(TE_sig, 3);
count  = sum(mpv1b < alpha, 3);

TE_avg = TE_sum ./ (count + eps);

%% =========================
% 3. REMOVE NEGATIVE VALUES
%% =========================

TE_avg(TE_avg < 0) = 0;

%% =========================
% 4. DIRECTIONALITY MATRIX
%% =========================

TE_dir = TE_avg - TE_avg';

%% =========================
% 5. VISUALIZATION
%% =========================

figure;
subplot(1,2,1);
imagesc(TE_avg);
colorbar;
title('Avg TE (Significant)');
axis square;

subplot(1,2,2);
imagesc(TE_dir);
colorbar;
title('Directional TE (i→j minus j→i)');
axis square;