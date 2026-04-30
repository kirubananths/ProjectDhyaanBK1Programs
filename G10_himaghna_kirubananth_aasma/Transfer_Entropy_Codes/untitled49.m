protocols = {'G1','G2'};
elecList = {'Oz', 'O1','O2','Fz','Pz'};
Fs = 1000;

Subject = '013AR';

path_file_G1 = 'D:\Users\JHS\Sem2\NSP\Total Dhyaan BK Data\Aasma\013AR\EEG\280122\G1\segmentedData\LFP\';
path_file_G2 = 'D:\Users\JHS\Sem2\NSP\Total Dhyaan BK Data\Aasma\013AR\EEG\280122\G2\segmentedData\LFP\';


data_G2_bl = cell(1,length(elecList));
data_G2_st = cell(1,length(elecList));

data_G2_gamma_bl = cell(1,length(elecList));
data_G2_gamma_st = cell(1,length(elecList));

data_G1_bl = cell(1,length(elecList));
data_G1_st = cell(1,length(elecList));

data_G1_gamma_bl = cell(1,length(elecList));
data_G1_gamma_st = cell(1,length(elecList));

for i=1:64
    filename = fullfile(path_file_G1,"elec"+i+".mat");
    loadedData = load(filename);
    if ismember(loadedData.analogInfo.labels, elecList)

        data_G1_bl{i} = normalize(loadedData.analogData(:,1:1250)); 
        data_G1_st{i} = normalize(loadedData.analogData(:,1251:2500)); 

        data_G1_gamma = normalize(bandpass(loadedData.analogData,[30 80],Fs));

        data_G1_gamma_bl{i} = data_G1_gamma(:,1:1250);
        data_G1_gamma_st{i} = data_G1_gamma(:,1251:2500);
    end
end

for i=1:64
    filename = fullfile(path_file_G2,"elec"+i+".mat");
    loadedData = load(filename);
    if ismember(loadedData.analogInfo.labels, elecList)

        data_G2_bl{i} = normalize(loadedData.analogData(:,1:1250)); 
        data_G2_st{i} = normalize(loadedData.analogData(:,1251:2500)); 

        data_G2_gamma = normalize(bandpass(loadedData.analogData,[30 80],Fs));

        data_G2_gamma_bl{i} = data_G2_gamma(:,1:1250);
        data_G2_gamma_st{i} = data_G2_gamma(:,1251:2500);
    end
end
disp('Done')

elec_indices = [2 13 16 17 18];

%parameters
params.nBins = 10;
params.k = 3;
params.l = 3;
params.tau = 2;

kNN = 4;
nSurr = 50;
trial = 120;
nelec = length(elec_indices);
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

for i=1:length(elec_indices)
    for j=1:length(elec_indices)
      
            
        for k=1:trial
            
            X1b = data_G1_bl{elec_indices(i)}(k,:)';
            Y1b = data_G1_bl{elec_indices(j)}(k,:)';

            X1s = data_G1_st{elec_indices(i)}(k,:)';
            Y1s = data_G1_st{elec_indices(j)}(k,:)';

            X2b = data_G2_bl{elec_indices(i)}(k,:)'; 
            Y2b = data_G2_bl{elec_indices(j)}(k,:)';

            X2s = data_G2_st{elec_indices(i)}(k,:)';
            Y2s = data_G2_st{elec_indices(j)}(k,:)';

            X1gb = data_G1_gamma_bl{elec_indices(i)}(k,:)';
            Y1gb = data_G1_gamma_bl{elec_indices(j)}(k,:)';

            X1gs = data_G1_gamma_st{elec_indices(i)}(k,:)'; 
            Y1gs = data_G1_gamma_st{elec_indices(j)}(k,:)';

            X2gb = data_G2_gamma_bl{elec_indices(i)}(k,:)'; 
            Y2gb = data_G2_gamma_bl{elec_indices(j)}(k,:)';

            X2gs = data_G2_gamma_st{elec_indices(i)}(k,:)'; 
            Y2gs = data_G2_gamma_st{elec_indices(j)}(k,:)';

            [TE1b, Co1b, pval1b] = discrete_pipeline(X1b, Y1b, params, nSurr);
            [TE1s, Co1s, pval1s] = discrete_pipeline(X1s, Y1s, params, nSurr);

            [TE2b, Co2b, pval2b] = discrete_pipeline(X2b, Y2b, params, nSurr);
            [TE2s, Co2s, pval2s] = discrete_pipeline(X2s, Y2s, params, nSurr);

            [TE1gs, Co1gs, pval1gs] = discrete_pipeline(X1gs, Y1gs, params, nSurr);
            [TE1gb, Co1gb, pval1gb] = discrete_pipeline(X1gb, Y1gb, params, nSurr);

            [TE2gb, Co2gb, pval2gb] = discrete_pipeline(X2gb, Y2gb, params, nSurr);
            [TE2gs, Co2gs, pval2gs] = discrete_pipeline(X2gs, Y2gs, params, nSurr);

            mTE1b(i,j,k) = TE1b; mTE1s(i,j,k) = TE1s;
            mTE2b(i,j,k) = TE2b; mTE2s(i,j,k) = TE2s;
            
            mTE1gb(i,j,k) = TE1gb; mTE1gs(i,j,k) = TE1gs;
            mTE2gb(i,j,k) = TE2gb; mTE2gs(i,j,k) = TE2gs;

            mCo1b(i,j,k) = Co1b; mCo1s(i,j,k) = Co1s;
            mCo2b(i,j,k) = Co2b; mCo2s(i,j,k) = Co2s;

            mCo1gb(i,j,k) = Co1gb; mCo1gs(i,j,k) = Co1gs;
            mCo2gb(i,j,k) = Co2gb; mCo2gs(i,j,k) = Co2gs;

            mpv1b(i,j,k) = pval1b; mpv1s(i,j,k) = pval1s;
            mpv2b(i,j,k) = pval2b; mpv2s(i,j,k) = pval2s;

            mpv1gb(i,j,k) = pval1gb; mpv1gs(i,j,k) = pval1gs;
            mpv2gb(i,j,k) = pval2gb; mpv2gs(i,j,k) = pval2gs;

        end
    end
end