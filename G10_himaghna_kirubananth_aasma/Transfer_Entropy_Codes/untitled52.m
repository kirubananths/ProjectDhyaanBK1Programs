

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
% 5. ENHANCED VISUALIZATION
%% =========================

labels = {'Fz','Pz','O1','Oz','O2'}; % match your electrode order

figure;

% ---------- (A) TE_avg ----------
subplot(1,3,1);
imagesc(TE_avg);
colorbar;

title({'Avg Significant TE', ...
       'Only edges with p < 0.05 retained', ...
       'Magnitude = information flow strength'});

axis square;

set(gca, 'XTick', 1:nElec, 'XTickLabel', labels, ...
         'YTick', 1:nElec, 'YTickLabel', labels);

xlabel('Target (j)');
ylabel('Source (i)');

% diagonal should be zero
hold on;
plot(1:nElec,1:nElec,'kx','MarkerSize',8,'LineWidth',1.5);
hold off;


% ---------- (B) Directionality ----------
subplot(1,3,2);
imagesc(TE_dir);
colorbar;

title({'Directional TE', ...
       'TE(i→j) - TE(j→i)', ...
       'Positive: i drives j | Negative: j drives i'});

axis square;

set(gca, 'XTick', 1:nElec, 'XTickLabel', labels, ...
         'YTick', 1:nElec, 'YTickLabel', labels);

xlabel('Target (j)');
ylabel('Source (i)');

colormap(gca, 'jet'); % diverging would be better, but MATLAB default ok


% ---------- (C) SIGNIFICANCE COUNT ----------
subplot(1,3,3);

sig_count = sum(mpv1b < alpha, 3);

imagesc(sig_count);
colorbar;

title({'Consistency across trials', ...
       '# trials where TE is significant'});

axis square;

set(gca, 'XTick', 1:nElec, 'XTickLabel', labels, ...
         'YTick', 1:nElec, 'YTickLabel', labels);

xlabel('Target (j)');
ylabel('Source (i)');