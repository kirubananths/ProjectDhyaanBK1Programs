protocols = {'G1','G2'};
elecList = {'Oz', 'O1','O2','Fz','Pz'};
Fs = 1000;

Subject = '013AR';

path_file_G1 = 'D:\Users\JHS\Sem2\NSP\Total Dhyaan BK Data\Aasma\013AR\EEG\280122\G1\segmentedData\LFP\';
path_file_G2 = 'D:\Users\JHS\Sem2\NSP\Total Dhyaan BK Data\Aasma\013AR\EEG\280122\G2\segmentedData\LFP\';

elec_indices = [2 13 16 17 18];
nelec = length(elec_indices);

data_G1_bl   = cell(1, nelec);  data_G1_st   = cell(1, nelec);
data_G1_gamma_bl = cell(1, nelec); data_G1_gamma_st = cell(1, nelec);
data_G2_bl   = cell(1, nelec);  data_G2_st   = cell(1, nelec);
data_G2_gamma_bl = cell(1, nelec); data_G2_gamma_st = cell(1, nelec);


for ei = 1:nelec
    i = elec_indices(ei);
    for gi = 1:2
        if gi == 1; path_file = path_file_G1; else; path_file = path_file_G2; end
        filename = fullfile(path_file, "elec" + i + ".mat");
        loadedData = load(filename);
        
        normData  = normalize(loadedData.analogData);
        gammaData = normalize(bandpass(loadedData.analogData, [30 80], Fs));
        
        if gi == 1
            data_G1_bl{ei}       = normData(:, 1:1250);
            data_G1_st{ei}       = normData(:, 1251:2500);
            data_G1_gamma_bl{ei} = gammaData(:, 1:1250);
            data_G1_gamma_st{ei} = gammaData(:, 1251:2500);
        else
            data_G2_bl{ei}       = normData(:, 1:1250);
            data_G2_st{ei}       = normData(:, 1251:2500);
            data_G2_gamma_bl{ei} = gammaData(:, 1:1250);
            data_G2_gamma_st{ei} = gammaData(:, 1251:2500);
        end
    end
end
disp('Done')

%elec_indices = [2 13 16 17 18];

%parameters
params.nBins = 8;
params.k = 3;
params.l = 3;
params.tau = 2;

kNN = 4;
nSurr = 50;
trial = 30;
%nelec = length(elec_indices);
nelpairs = nelec * (nelec-1);

mTE1b = zeros(nelec,nelec,trial); mCo1b = zeros(nelec,nelec,trial); mpv1b = zeros(nelec,nelec,trial);
mTE1s = zeros(nelec,nelec,trial); mCo1s = zeros(nelec,nelec,trial); mpv1s = zeros(nelec,nelec,trial);

mTE1gb = zeros(nelec,nelec,trial); mCo1gb = zeros(nelec,nelec,trial); mpv1gb = zeros(nelec,nelec,trial);
mTE1gs = zeros(nelec,nelec,trial); mCo1gs = zeros(nelec,nelec,trial); mpv1gs = zeros(nelec,nelec,trial);

mTE2b = zeros(nelec,nelec,trial); mCo2b = zeros(nelec,nelec,trial); mpv2b = zeros(nelec,nelec,trial);
mTE2s = zeros(nelec,nelec,trial); mCo2s = zeros(nelec,nelec,trial); mpv2s = zeros(nelec,nelec,trial);

mTE2gb = zeros(nelec,nelec,trial); mCo2gb = zeros(nelec,nelec,trial); mpv2gb = zeros(nelec,nelec,trial);
mTE2gs = zeros(nelec,nelec,trial); mCo2gs = zeros(nelec,nelec,trial); mpv2gs = zeros(nelec,nelec,trial);


TE_corr = zeros(nelec,nelec,trial);
pval = zeros(nelec,nelec,trial);

for i=1:nelec
    for j=1:nelec

        fprintf('Processing pair (%d,%d) of (%d,%d)\n', i, j, nelec, nelec);

        if i==j; continue; end

        Xi1b = data_G1_bl{i};    Xj1b = data_G1_bl{j};
        Xi1s = data_G1_st{i};    Xj1s = data_G1_st{j};
        Xi2b = data_G2_bl{i};    Xj2b = data_G2_bl{j};
        Xi2s = data_G2_st{i};    Xj2s = data_G2_st{j};
        Xi1gb = data_G1_gamma_bl{i}; Xj1gb = data_G1_gamma_bl{j};
        Xi1gs = data_G1_gamma_st{i}; Xj1gs = data_G1_gamma_st{j};
        Xi2gb = data_G2_gamma_bl{i}; Xj2gb = data_G2_gamma_bl{j};
        Xi2gs = data_G2_gamma_st{i}; Xj2gs = data_G2_gamma_st{j};

        % Temp arrays for parfor (can't index into 3D arrays directly)
        te1b=zeros(1,trial); co1b=zeros(1,trial); pv1b=zeros(1,trial);
        te1s=zeros(1,trial); co1s=zeros(1,trial); pv1s=zeros(1,trial);
        te2b=zeros(1,trial); co2b=zeros(1,trial); pv2b=zeros(1,trial);
        te2s=zeros(1,trial); co2s=zeros(1,trial); pv2s=zeros(1,trial);
        te1gb=zeros(1,trial); co1gb=zeros(1,trial); pv1gb=zeros(1,trial);
        te1gs=zeros(1,trial); co1gs=zeros(1,trial); pv1gs=zeros(1,trial);
        te2gb=zeros(1,trial); co2gb=zeros(1,trial); pv2gb=zeros(1,trial);
        te2gs=zeros(1,trial); co2gs=zeros(1,trial); pv2gs=zeros(1,trial);
            
        parfor k = 1:trial
            [t, c, p] = discrete_pipeline(Xi1b(k,:)', Xj1b(k,:)', params, nSurr);
    te1b(k) = t; co1b(k) = c; pv1b(k) = p;

    [t, c, p] = discrete_pipeline(Xi1s(k,:)', Xj1s(k,:)', params, nSurr);
    te1s(k) = t; co1s(k) = c; pv1s(k) = p;

    [t, c, p] = discrete_pipeline(Xi2b(k,:)', Xj2b(k,:)', params, nSurr);
    te2b(k) = t; co2b(k) = c; pv2b(k) = p;

    [t, c, p] = discrete_pipeline(Xi2s(k,:)', Xj2s(k,:)', params, nSurr);
    te2s(k) = t; co2s(k) = c; pv2s(k) = p;

    [t, c, p] = discrete_pipeline(Xi1gb(k,:)', Xj1gb(k,:)', params, nSurr);
    te1gb(k) = t; co1gb(k) = c; pv1gb(k) = p;

    [t, c, p] = discrete_pipeline(Xi1gs(k,:)', Xj1gs(k,:)', params, nSurr);
    te1gs(k) = t; co1gs(k) = c; pv1gs(k) = p;

    [t, c, p] = discrete_pipeline(Xi2gb(k,:)', Xj2gb(k,:)', params, nSurr);
    te2gb(k) = t; co2gb(k) = c; pv2gb(k) = p;

    [t, c, p] = discrete_pipeline(Xi2gs(k,:)', Xj2gs(k,:)', params, nSurr);
    te2gs(k) = t; co2gs(k) = c; pv2gs(k) = p;
        end
            mTE1b(i,j,:)=te1b; mCo1b(i,j,:)=co1b; mpv1b(i,j,:)=pv1b;
            mTE1s(i,j,:)=te1s; mCo1s(i,j,:)=co1s; mpv1s(i,j,:)=pv1s;
            mTE2b(i,j,:)=te2b; mCo2b(i,j,:)=co2b; mpv2b(i,j,:)=pv2b;
            mTE2s(i,j,:)=te2s; mCo2s(i,j,:)=co2s; mpv2s(i,j,:)=pv2s;
            mTE1gb(i,j,:)=te1gb; mCo1gb(i,j,:)=co1gb; mpv1gb(i,j,:)=pv1gb;
            mTE1gs(i,j,:)=te1gs; mCo1gs(i,j,:)=co1gs; mpv1gs(i,j,:)=pv1gs;
            mTE2gb(i,j,:)=te2gb; mCo2gb(i,j,:)=co2gb; mpv2gb(i,j,:)=pv2gb;
            mTE2gs(i,j,:)=te2gs; mCo2gs(i,j,:)=co2gs; mpv2gs(i,j,:)=pv2gs;
     end
end