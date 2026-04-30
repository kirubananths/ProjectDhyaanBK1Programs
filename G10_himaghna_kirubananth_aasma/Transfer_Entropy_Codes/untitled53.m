labels = {'Fz','Pz','O1','Oz','O2'};
alpha = 0.05;

function TE_avg = compute_TE_avg(mCo, mpv, alpha)

TE_sig = mCo .* (mpv < alpha);
TE_sum = sum(TE_sig, 3);
count  = sum(mpv < alpha, 3);

TE_avg = TE_sum ./ (count + eps);
TE_avg(TE_avg < 0) = 0;

end


TE_G1 = compute_TE_avg(mCo1s, mpv1s, alpha);
TE_G2 = compute_TE_avg(mCo2s, mpv2s, alpha);

TE_diff = TE_G2 - TE_G1;


figure;

subplot(1,3,1);
imagesc(TE_G1);
colorbar;
title('G1 (Meditation)');
axis square;

subplot(1,3,2);
imagesc(TE_G2);
colorbar;
title('G2 (Control)');
axis square;

subplot(1,3,3);
imagesc(TE_diff);
colorbar;
title('G2 - G1 (Meditation Effect)');
axis square;

for k = 1:3
    subplot(1,3,k);
    set(gca,'XTick',1:5,'XTickLabel',labels,...
            'YTick',1:5,'YTickLabel',labels);
    xlabel('Target'); ylabel('Source');
end

TE_bl = compute_TE_avg(mCo1b, mpv1b, alpha);
TE_st = compute_TE_avg(mCo1s, mpv1s, alpha);

TE_diff = TE_st - TE_bl;

figure;

subplot(1,3,1);
imagesc(TE_bl);
colorbar;
title('Baseline (G1)');
axis square;

subplot(1,3,2);
imagesc(TE_st);
colorbar;
title('Stimulus (G1)');
axis square;

subplot(1,3,3);
imagesc(TE_diff);
colorbar;
title('Stim - Base (G1)');
axis square;


TE_bl2 = compute_TE_avg(mCo2b, mpv2b, alpha);
TE_st2 = compute_TE_avg(mCo2s, mpv2s, alpha);

TE_diff2 = TE_st2 - TE_bl2;

figure;

subplot(1,3,1);
imagesc(TE_bl2);
colorbar;
title('Baseline (G2)');
axis square;

subplot(1,3,2);
imagesc(TE_st2);
colorbar;
title('Stimulus (G2)');
axis square;

subplot(1,3,3);
imagesc(TE_diff2);
colorbar;
title('Stim - Base (G2)');
axis square;



TE_raw   = compute_TE_avg(mCo1s, mpv1s, alpha);
TE_gamma = compute_TE_avg(mCo1gs, mpv1gs, alpha);

TE_diff = TE_gamma - TE_raw;

figure;

subplot(1,3,1);
imagesc(TE_raw);
colorbar;
title('Broadband');
axis square;

subplot(1,3,2);
imagesc(TE_gamma);
colorbar;
title('Gamma (30–80 Hz)');
axis square;

subplot(1,3,3);
imagesc(TE_diff);
colorbar;
title('Gamma - Broadband');
axis square;

