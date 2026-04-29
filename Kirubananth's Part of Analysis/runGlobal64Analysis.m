function runGlobal64Analysis(varargin)
% RUNGLOBAL64ANALYSIS  Paired group-level analysis of meditator vs
%                      non-meditator EEG using all 64 channels.
% =========================================================================
%
%   Three core metrics are computed for every paired subject:
%
%       LZC : Lempel-Ziv Complexity (broadband signal complexity)
%       PAC : Phase-Amplitude Coupling - quantified by the Modulation
%             Index of Tort et al. 2010 (theta-gamma)
%       FC  : Functional Connectivity quantified by Mutual Information
%             (full Neeg x Neeg matrix, in bits)
%
%   The script uses the FULL 64-channel EEG montage at both per-pair and
%   group level - the bad-electrode flagging is reported for diagnostic
%   transparency but does NOT exclude any channel from the calculation.
%   Per-trial QC (badTrials) still applies.
%
%   USAGE
%   -----
%       cd /path/to/Meditators
%       runGlobal64Analysis;                              % default EC2
%       runGlobal64Analysis('protocol', 'EO1');
%       runGlobal64Analysis('refreshCache', true);
%       runGlobal64Analysis('makePerPairFigures', false); % skip 30 pair figs
%       runGlobal64Analysis('runFC', false);              % skip FC (~saves 1 hr)
%       runGlobal64Analysis('statsTest', 'wilcoxon');
%
%   FOLDER LAYOUT
%   -------------
%       <study root>/
%           Meditators/   (paste this script HERE)
%               <medID>/...
%               pairList.csv
%           Non-Meditators/
%               <nonID>/...
%
%   Supported per-subject layouts (probed in this order):
%       1. <id>/<protocol>/                             (FLAT)
%       2. <id>/<protocol>/LFP/
%       3. <id>/<protocol>/segmentedData/LFP/
%       4. <id>/EEG/<date>/<protocol>/segmentedData/LFP/
%
%   OUTPUTS in Results_Global64_<protocol>/
%   ---------------------------------------
%       cache/pair_<medID>_<nonID>.mat       per-pair raw results
%       GLOBAL64_RESULTS.mat                 everything aggregated
%       global64_summary_stats.csv           scalar metrics summary
%       global64_FC_pairs.csv                per-connection summary
%
%       Fig_LZC_topoplots.png/.fig           Med | Non | Paired Delta
%       Fig_LZC_summary.png/.fig             4-panel deep dive
%       Fig_PAC_topoplots.png/.fig           Med | Non | Paired Delta
%       Fig_PAC_summary.png/.fig             4-panel deep dive
%       Fig_FC_matrices.png/.fig             Mutual Information matrices
%       Fig_FC_summary.png/.fig              top-K connections + hubs + ...
%       Fig_FC_network.png/.fig              head-outline network
%       Fig_GlobalSummary.png/.fig           one-figure overview
%       pairs/pair_<medID>_<nonID>.png/.fig  one figure per pair
%
%   Author: Kirubananth S
% =========================================================================


%% ========================================================================
% CONFIG
% =========================================================================
cfg.medRoot          = pwd;
cfg.nonRoot          = fullfile(pwd, '..', 'Non-Meditators');
cfg.protocol         = 'EC2';
cfg.timeWindow       = [-1.0 1.0];
cfg.outputDir        = '';                % auto: pwd/Results_Global64_<proto>
cfg.cacheDir         = '';                % auto: outputDir/cache

cfg.refreshCache     = false;
cfg.reuseGroupCache  = true;
cfg.groupCacheDir    = '';                % auto: ../Results_Group_<proto>/cache

% Which analyses to run.
cfg.runLZC           = true;
cfg.runPAC           = true;
cfg.runFC            = true;              % Mutual-Information functional connectivity

% Algorithm parameters.
cfg.lzcBand          = [1 100];
cfg.thetaRange       = [4 8];
cfg.gammaRange       = [30 80];
cfg.nPhaseBins       = 18;
cfg.fcBand           = [1 45];            % bandpass before mutual information
cfg.fcBins           = 16;                % quantile bins for the MI estimator

% Statistics.
cfg.statsTest        = 'pairedT';         % 'pairedT' or 'wilcoxon'
cfg.alpha            = 0.05;

% Visualisation.
cfg.makePerPairFigures = true;
cfg.topKConnections  = 30;
cfg.topKElectrodes   = 12;

% Color limits for connectivity matrices: percentile-based so the diagonal
% (MI of an electrode with itself = NaN) and a few tail outliers don't
% dominate the colormap.  98th percentile is conservative.
cfg.fcMatrixPctile   = 98;

% Electrode selection.  Three modes are available:
%   useAllElectrodes = true            -> use all 64 channels (default)
%   useAllElectrodes = false           -> use channels clean in
%                                          >= minGoodFraction of subjects
%   useAllElectrodes = false,
%   minGoodFraction  = 1.0             -> strict intersection (clean in
%                                          every subject of every pair)
% In every mode the per-trial QC (badTrials) is still applied.
cfg.useAllElectrodes = true;              % Set this to 'false' if you want to use those electrodes that are good in atleast 80 percent of the Total Subjects (Meditators + Controls)
cfg.minGoodFraction  = 0.80;
cfg.minGlobalElecs   = 5;                 % warn if final selection < this

% Apply name-value overrides.
for k = 1:2:numel(varargin)
    cfg.(varargin{k}) = varargin{k+1};
end

% Resolve auto-fill paths.
if isempty(cfg.outputDir)
    cfg.outputDir = fullfile(pwd, sprintf('Results_Global64_%s', cfg.protocol));
end
if isempty(cfg.cacheDir)
    cfg.cacheDir = fullfile(cfg.outputDir, 'cache');
end
if isempty(cfg.groupCacheDir)
    cfg.groupCacheDir = fullfile(pwd, ...
        sprintf('Results_Group_%s', cfg.protocol), 'cache');
end
pairsDir = fullfile(cfg.outputDir, 'pairs');
mkdirSafe(cfg.outputDir);
mkdirSafe(cfg.cacheDir);
if cfg.makePerPairFigures, mkdirSafe(pairsDir); end


%% ========================================================================
% LOCATE eegUtils
% =========================================================================
if exist('eegUtils', 'file') ~= 2
    candidates = {pwd, fullfile(pwd,'..','Codes'), fullfile(pwd,'..'), ...
                  fileparts(mfilename('fullpath'))};
    for k = 1:numel(candidates)
        if exist(fullfile(candidates{k}, 'eegUtils.m'), 'file')
            addpath(candidates{k}); break;
        end
    end
end
if exist('eegUtils', 'file') ~= 2
    error('eegUtils.m not found.');
end
u = eegUtils();


%% ========================================================================
% PAIR LIST
% =========================================================================
csvPath = fullfile(fileparts(mfilename('fullpath')), 'pairList.csv');
if ~exist(csvPath, 'file')
    csvPath = fullfile(pwd, 'pairList.csv');
end
if ~exist(csvPath, 'file')
    error(['pairList.csv not found.  See pairList.csv in the parent ' ...
           'folder for format.']);
end
cfg.pairs = readPairList(csvPath);
keep = ~cellfun(@isempty, cfg.pairs(:,1)) & ...
       ~cellfun(@isempty, cfg.pairs(:,2));
cfg.pairs = cfg.pairs(keep, :);
nPairs = size(cfg.pairs, 1);
if nPairs == 0, error('No usable pairs.'); end


%% ========================================================================
% HEADER
% =========================================================================
u.printDivider('GLOBAL 64-CHANNEL PAIRED ANALYSIS (Med vs Non-Med)');
fprintf('Med root         : %s\n', cfg.medRoot);
fprintf('Non-Med root     : %s\n', cfg.nonRoot);
fprintf('Protocol         : %s\n', cfg.protocol);
fprintf('Time window      : [%.3f %.3f] s\n', cfg.timeWindow);
fprintf('Pairs            : %d\n', nPairs);
fprintf('Analyses         : LZC=%d  PAC=%d  FC=%d\n', cfg.runLZC, cfg.runPAC, cfg.runFC);
fprintf('Stats            : %s   alpha=%.3f (FDR-BH)\n', cfg.statsTest, cfg.alpha);
fprintf('Output           : %s\n', cfg.outputDir);


%% ========================================================================
% PER-PAIR PROCESSING
% =========================================================================
pairResults = cell(nPairs, 1);
nReused = 0;  nFresh = 0;  nOwnCache = 0;

for p = 1:nPairs
    [medID, nonID, medDate, nonDate] = unpackPairRow(cfg.pairs, p);
    fprintf('\n[Pair %2d/%2d]  Med:%s  Non:%s\n', p, nPairs, medID, nonID);

    [R, source] = loadOrComputePair(medID, nonID, medDate, nonDate, cfg, u);
    if isempty(R)
        fprintf('   SKIP - data not found\n');
        continue;
    end
    pairResults{p} = R;
    switch source
        case 'own',     nOwnCache = nOwnCache + 1; fprintf('   own-cache hit\n');
        case 'group',   nReused   = nReused + 1;   fprintf('   reused runGroupAnalysis cache\n');
        case 'fresh',   nFresh    = nFresh + 1;
    end
end

ok = ~cellfun(@isempty, pairResults);
pairResults = pairResults(ok);
nPairs = numel(pairResults);
if nPairs == 0, error('No pairs loaded successfully.'); end

fprintf('\n%d pair(s) ready  (own-cache=%d, group-cache=%d, fresh=%d)\n', ...
        nPairs, nOwnCache, nReused, nFresh);


%% ========================================================================
% ELECTRODE SET DIAGNOSTIC
% =========================================================================
% All 64 electrodes are used, but we still print which ones are flagged
% bad and which subjects are outliers - useful for interpreting results.
Neeg = numel(pairResults{1}.eegLabels);
nSubjects = 2 * nPairs;
labels = pairResults{1}.eegLabels;

electrodeGoodness = zeros(Neeg, 1);
for p = 1:nPairs
    electrodeGoodness(pairResults{p}.medGoodElecs) = ...
        electrodeGoodness(pairResults{p}.medGoodElecs) + 1;
    electrodeGoodness(pairResults{p}.nonGoodElecs) = ...
        electrodeGoodness(pairResults{p}.nonGoodElecs) + 1;
end

% Apply electrode-selection rule.
if cfg.useAllElectrodes
    globalElecs = (1:Neeg)';
    selectionMode = 'all 64 channels (per-electrode QC bypassed)';
else
    threshold   = ceil(cfg.minGoodFraction * nSubjects);
    globalElecs = find(electrodeGoodness >= threshold);
    selectionMode = sprintf('clean in >= %.0f%% of subjects (>= %d / %d)', ...
                            100*cfg.minGoodFraction, threshold, nSubjects);
end
nG = numel(globalElecs);

fprintf('\n=== Electrode selection ===\n');
fprintf('Criterion       : %s\n', selectionMode);
fprintf('Selected        : %d / %d channels  (%.1f%%)\n', ...
        nG, Neeg, 100*nG/Neeg);
if ~cfg.useAllElectrodes && nG > 0
    fprintf('Indices         : %s\n', mat2str(globalElecs(:)'));
    fprintf('Labels          : %s\n', strjoin(labels(globalElecs), ', '));
end
fprintf('Per-trial QC (badTrials) is applied in every mode.\n');

[gSorted, gIdx] = sort(electrodeGoodness, 'descend');
fprintf('\nMost-reliable channels by per-subject QC (top 10):\n');
for k = 1:min(10, Neeg)
    fprintf('  %-4s (idx %2d): clean in %2d / %d subjects (%3.0f%%)\n', ...
            labels{gIdx(k)}, gIdx(k), gSorted(k), nSubjects, ...
            100*gSorted(k)/nSubjects);
end
fprintf('Least-reliable channels by per-subject QC (bottom 10):\n');
for k = max(1, Neeg-9):Neeg
    fprintf('  %-4s (idx %2d): clean in %2d / %d subjects (%3.0f%%)\n', ...
            labels{gIdx(k)}, gIdx(k), gSorted(k), nSubjects, ...
            100*gSorted(k)/nSubjects);
end

subjBad = zeros(nSubjects, 1);
subjIDs = cell(nSubjects, 1);
for p = 1:nPairs
    subjBad(2*p-1) = Neeg - numel(pairResults{p}.medGoodElecs);
    subjIDs{2*p-1} = ['M:' pairResults{p}.medID];
    subjBad(2*p)   = Neeg - numel(pairResults{p}.nonGoodElecs);
    subjIDs{2*p}   = ['N:' pairResults{p}.nonID];
end
[bSorted, bIdx] = sort(subjBad, 'descend');
fprintf('\nSubjects with most QC-flagged electrodes (top 10):\n');
for k = 1:min(10, nSubjects)
    fprintf('  %-12s : %2d flagged channels\n', subjIDs{bIdx(k)}, bSorted(k));
end

diagInfo.electrodeGoodness = electrodeGoodness;
diagInfo.subjBadCount      = subjBad;
diagInfo.subjIDs           = subjIDs;
diagInfo.nSubjects         = nSubjects;


%% ========================================================================
% AGGREGATION
% =========================================================================
agg = struct();
agg.cfg          = cfg;
agg.refLabels    = pairResults{1}.eegLabels;
agg.refXY        = pairResults{1}.elecXY;
agg.globalElecs  = globalElecs;
agg.nPairs       = nPairs;
agg.windowLabel  = pairResults{1}.windowLabel;
agg.diagnostic   = diagInfo;
agg.pairIDs      = cellfun(@(r) sprintf('%s_%s', r.medID, r.nonID), ...
                           pairResults, 'UniformOutput', false);

if cfg.runLZC
    agg.LZC = aggregateScalar(pairResults, 'medLZC', 'nonLZC', globalElecs, cfg);
end
if cfg.runPAC
    agg.PAC = aggregateScalar(pairResults, 'medPAC', 'nonPAC', globalElecs, cfg);
end
if cfg.runFC
    agg.FC  = aggregateMatrix(pairResults, 'medFC', 'nonFC', globalElecs, cfg);
end

save(fullfile(cfg.outputDir, 'GLOBAL64_RESULTS.mat'), ...
     'cfg', 'agg', 'pairResults', '-v7.3');
fprintf('\nSaved: %s\n', fullfile(cfg.outputDir, 'GLOBAL64_RESULTS.mat'));


%% ========================================================================
% CSV SUMMARIES
% =========================================================================
writeScalarCSV(agg, cfg);
if cfg.runFC
    writeMatrixCSV(agg, cfg);
end


%% ========================================================================
% PLOTS
% =========================================================================
fprintf('\n--- Generating figures ---\n');

if cfg.runLZC
    plotMetricTopoplots(agg, 'LZC', 'Lempel-Ziv Complexity', ...
                        'normalised LZC', cfg, u);
    plotMetricSummary  (agg, 'LZC', 'Lempel-Ziv Complexity', ...
                        'normalised LZC', cfg, u);
end
if cfg.runPAC
    pacName = sprintf('Theta(%g-%g) - Gamma(%g-%g) PAC', ...
                       cfg.thetaRange(1), cfg.thetaRange(2), ...
                       cfg.gammaRange(1), cfg.gammaRange(2));
    plotMetricTopoplots(agg, 'PAC', pacName, 'Modulation Index', cfg, u);
    plotMetricSummary  (agg, 'PAC', pacName, 'Modulation Index', cfg, u);
end
if cfg.runFC
    plotFCMatrices(agg, cfg, u);
    plotFCSummary (agg, cfg, u);
    plotFCNetwork (agg, cfg, u);
end

plotGlobalSummary(agg, cfg, u);


%% ========================================================================
% PER-PAIR FIGURES
% =========================================================================
if cfg.makePerPairFigures
    fprintf('\n--- Per-pair figures (%d pairs) ---\n', nPairs);
    for p = 1:nPairs
        plotPairFigure(pairResults{p}, agg, p, nPairs, pairsDir, cfg, u);
    end
end

u.printDivider('GLOBAL 64-CHANNEL ANALYSIS COMPLETE');
fprintf('Results in: %s\n', cfg.outputDir);
fprintf('All figures left open on screen for manual export.\n');
end


%% ========================================================================
% =================  PER-PAIR PROCESSING  ================================
%% ========================================================================

function [R, source] = loadOrComputePair(medID, nonID, medDate, nonDate, cfg, u)
% Try own cache, then group-cache (legacy), then compute fresh.
R = []; source = '';
fname = sprintf('pair_%s_%s.mat', medID, nonID);

% --- 1. Own cache ----------------------------------------------------
ownPath = fullfile(cfg.cacheDir, fname);
if exist(ownPath, 'file') && ~cfg.refreshCache
    L = load(ownPath, 'R');
    if pairCacheIsValid(L.R, cfg)
        R = L.R; source = 'own'; return;
    end
end

% --- 2. Legacy cache from runGroupAnalysis ---------------------------
% Only useful when the legacy file used the same field names.  Modern
% caches use medFC/nonFC; legacy ones used medMI/nonMI.  We rename
% on the fly when migrating.
if cfg.reuseGroupCache && ~cfg.refreshCache
    altPath = fullfile(cfg.groupCacheDir, fname);
    if exist(altPath, 'file')
        L = load(altPath, 'R');
        Rmig = migrateLegacyFields(L.R);
        if pairCacheIsValid(Rmig, cfg)
            R = Rmig; source = 'group';
            save(ownPath, 'R'); %#ok<NASGU>
            return;
        else
            fprintf('   group cache present but missing required fields - recomputing\n');
        end
    end
end

% --- 3. Compute fresh ------------------------------------------------
medLFP = findLFPPath(cfg.medRoot, medID, cfg.protocol, medDate);
nonLFP = findLFPPath(cfg.nonRoot, nonID, cfg.protocol, nonDate);
if isempty(medLFP) || isempty(nonLFP), return; end
fprintf('   computing fresh...\n');
fprintf('     med LFP: %s\n', medLFP);
fprintf('     non LFP: %s\n', nonLFP);

Sm = u.loadSubject(medLFP);
Sn = u.loadSubject(nonLFP);

R = struct();
R.medID = medID;        R.nonID = nonID;
R.medFolder = medLFP;   R.nonFolder = nonLFP;
R.eegLabels = Sm.eegLabels;
R.elecXY    = Sm.elecXY;
R.medGoodElecs  = setdiff(find(Sm.elecAvailable), Sm.badElectrodes);
R.nonGoodElecs  = setdiff(find(Sn.elecAvailable), Sn.badElectrodes);
R.medHighPri    = Sm.highPriorityElecs;
R.nonHighPri    = Sn.highPriorityElecs;
R.bothGoodElecs = intersect(R.medGoodElecs, R.nonGoodElecs);

if cfg.runLZC
    fprintf('     LZC...\n');
    R.medLZC = computeLZC(Sm, cfg, u);
    R.nonLZC = computeLZC(Sn, cfg, u);
end
if cfg.runPAC
    fprintf('     PAC (Modulation Index)...\n');
    R.medPAC = computePAC(Sm, cfg, u);
    R.nonPAC = computePAC(Sn, cfg, u);
end
if cfg.runFC
    fprintf('     FC (Mutual Information)...\n');
    R.medFC = computeFC(Sm, cfg, u);
    R.nonFC = computeFC(Sn, cfg, u);
end

[~, ~, tInfoM] = u.getTimeWindowIdx(Sm.timeVals, cfg.timeWindow);
R.windowLabel = tInfoM.label;

save(ownPath, 'R');
source = 'fresh';
end


function Rout = migrateLegacyFields(R)
% Rename legacy medMI/nonMI -> medFC/nonFC (mutual-information matrices).
Rout = R;
if isfield(R, 'medMI') && ~isfield(R, 'medFC'), Rout.medFC = R.medMI; end
if isfield(R, 'nonMI') && ~isfield(R, 'nonFC'), Rout.nonFC = R.nonMI; end
end


function ok = pairCacheIsValid(R, cfg)
if ~isstruct(R), ok = false; return; end
required = {'bothGoodElecs', 'eegLabels', 'medGoodElecs', 'nonGoodElecs'};
for k = 1:numel(required)
    if ~isfield(R, required{k}), ok = false; return; end
end
if cfg.runLZC && (~isfield(R, 'medLZC') || ~isfield(R, 'nonLZC'))
    ok = false; return;
end
if cfg.runPAC && (~isfield(R, 'medPAC') || ~isfield(R, 'nonPAC'))
    ok = false; return;
end
if cfg.runFC && (~isfield(R, 'medFC') || ~isfield(R, 'nonFC'))
    ok = false; return;
end
ok = true;
end


%% ========================================================================
% =================  PER-SUBJECT COMPUTE HELPERS  ========================
%% ========================================================================

function lzc = computeLZC(S, cfg, u)
[tIdx, tWin, ~] = u.getTimeWindowIdx(S.timeVals, cfg.timeWindow);
nSamp = numel(tWin);
keepTrials = setdiff(1:S.totalTrials, S.badTrials(:)');
goodElecs  = pickGoodElecs(S, cfg);

Neeg = numel(S.eegLabels);
lzc  = NaN(Neeg, 1);
for e = goodElecs
    raw = double(squeeze(S.data(e, keepTrials, tIdx)));
    if isempty(raw), continue; end
    sig = u.bandpassFilter(raw', S.Fs, cfg.lzcBand(1), cfg.lzcBand(2))';
    sig = (sig - mean(sig,2)) ./ (std(sig,0,2) + eps);
    perTrial = NaN(size(sig,1), 1);
    for ti = 1:size(sig,1)
        b = double(sig(ti,:) > median(sig(ti,:)));
        c = lz76local(b);
        perTrial(ti) = c * log2(nSamp) / nSamp;
    end
    lzc(e) = mean(perTrial, 'omitnan');
end
end


function pac = computePAC(S, cfg, u)
% Tort 2010 Modulation Index: theta-phase / gamma-amplitude.
[tIdx, ~, ~] = u.getTimeWindowIdx(S.timeVals, cfg.timeWindow);
keepTrials = setdiff(1:S.totalTrials, S.badTrials(:)');
goodElecs  = pickGoodElecs(S, cfg);
Neeg = numel(S.eegLabels);
pac  = NaN(Neeg, 1);
for e = goodElecs
    sigMat = double(squeeze(S.data(e, keepTrials, :)));
    if isvector(sigMat), sigMat = sigMat(:)'; end
    thetaSig = u.bandpassFilter(sigMat', S.Fs, cfg.thetaRange(1), cfg.thetaRange(2))';
    gammaSig = u.bandpassFilter(sigMat', S.Fs, cfg.gammaRange(1), cfg.gammaRange(2))';
    [phs, amp] = trialwiseHilbertLocal(thetaSig, gammaSig, tIdx);
    pac(e) = tortMIlocal(phs, amp, cfg.nPhaseBins);
end
end


function FC = computeFC(S, cfg, u)
% Mutual-Information functional connectivity (bits).
[tIdx, tWin, ~] = u.getTimeWindowIdx(S.timeVals, cfg.timeWindow);
nWin = numel(tWin);
keepTrials = setdiff(1:S.totalTrials, S.badTrials(:)');
goodElecs  = pickGoodElecs(S, cfg);
nGood = numel(goodElecs);

% tIdx is a logical mask (length = numel(timeVals)), so use nWin (count
% of true entries = numel(tWin)) when allocating the per-electrode
% concatenated-signal matrix.
sigCat = zeros(nGood, nWin * numel(keepTrials));
for ei = 1:nGood
    e = goodElecs(ei);
    sigMat = double(squeeze(S.data(e, keepTrials, :)));
    if isvector(sigMat), sigMat = sigMat(:)'; end
    flt = bandpassRowsLocal(sigMat, S.Fs, cfg.fcBand(1), cfg.fcBand(2));
    flt = flt(:, tIdx);
    sigCat(ei,:) = flt(:)';
end

binIdx = zeros(size(sigCat),'uint8');
for ei = 1:nGood
    binIdx(ei,:) = quantileBinLocal(sigCat(ei,:), cfg.fcBins);
end

Neeg = numel(S.eegLabels);
FC = NaN(Neeg, Neeg);
for ai = 1:nGood-1
    for bi = ai+1:nGood
        m = miFromBinsLocal(binIdx(ai,:), binIdx(bi,:), cfg.fcBins);
        FC(goodElecs(ai), goodElecs(bi)) = m;
        FC(goodElecs(bi), goodElecs(ai)) = m;
    end
end
end


%% ========================================================================
% =================  AGGREGATION  ========================================
%% ========================================================================

function A = aggregateScalar(pairResults, fMed, fNon, gE, cfg)
nP = numel(pairResults);
nG = numel(gE);

medAll  = NaN(nP, nG);
nonAll  = NaN(nP, nG);
for p = 1:nP
    R = pairResults{p};
    medAll(p,:) = R.(fMed)(gE).';
    nonAll(p,:) = R.(fNon)(gE).';
end
diffAll = medAll - nonAll;

A.electrodes  = gE;
A.medAll      = medAll;
A.nonAll      = nonAll;
A.diffAll     = diffAll;
A.meanMed     = mean(medAll, 1, 'omitnan');
A.meanNon     = mean(nonAll, 1, 'omitnan');
A.meanDiff    = mean(diffAll, 1, 'omitnan');
A.semDiff     = std (diffAll, 0, 1, 'omitnan') ./ sqrt(max(1, sum(~isnan(diffAll), 1)));
A.cohenD      = A.meanDiff ./ (std(diffAll, 0, 1, 'omitnan') + eps);
A.tStat       = A.meanDiff ./ (A.semDiff + eps);

pVals = NaN(1, nG);
for k = 1:nG
    pVals(k) = pairedPvalue(diffAll(:, k), cfg.statsTest);
end
A.pVals    = pVals;
A.pAdj     = fdrBHlocal(pVals);
A.signMask = A.pAdj < cfg.alpha;
end


function A = aggregateMatrix(pairResults, fMed, fNon, gE, cfg)
nP = numel(pairResults);
nG = numel(gE);

medAll  = NaN(nP, nG, nG);
nonAll  = NaN(nP, nG, nG);
for p = 1:nP
    R = pairResults{p};
    medAll(p,:,:) = R.(fMed)(gE, gE);
    nonAll(p,:,:) = R.(fNon)(gE, gE);
end
diffAll = medAll - nonAll;

A.electrodes = gE;
A.medAll     = medAll; %#ok<STRNU>
A.nonAll     = nonAll;
A.diffAll    = diffAll;
A.meanMed    = squeeze(mean(medAll,  1, 'omitnan'));
A.meanNon    = squeeze(mean(nonAll,  1, 'omitnan'));
A.meanDiff   = squeeze(mean(diffAll, 1, 'omitnan'));
sd           = squeeze(std (diffAll, 0, 1, 'omitnan'));
A.semDiff    = sd ./ sqrt(max(1, squeeze(sum(~isnan(diffAll), 1))));
A.cohenD     = A.meanDiff ./ (sd + eps);

pVals = NaN(nG, nG);
for i = 1:nG-1
    for j = i+1:nG
        d = squeeze(diffAll(:, i, j));
        pp = pairedPvalue(d, cfg.statsTest);
        pVals(i, j) = pp;
        pVals(j, i) = pp;
    end
end
A.pVals = pVals;

upper = triu(true(nG), 1);
pVec  = pVals(upper);
pAdj  = fdrBHlocal(pVec);
pAdjMat = NaN(nG, nG);
pAdjMat(upper) = pAdj;
pAdjMat = pAdjMat + pAdjMat';
A.pAdj     = pAdjMat;
A.signMask = pAdjMat < cfg.alpha;

strMed = squeeze(sum(medAll, 3, 'omitnan'));
strNon = squeeze(sum(nonAll, 3, 'omitnan'));
strDiff = strMed - strNon;
A.nodeStrength.medAll   = strMed;
A.nodeStrength.nonAll   = strNon;
A.nodeStrength.diffAll  = strDiff;
A.nodeStrength.meanMed  = mean(strMed, 1, 'omitnan');
A.nodeStrength.meanNon  = mean(strNon, 1, 'omitnan');
A.nodeStrength.meanDiff = mean(strDiff, 1, 'omitnan');
A.nodeStrength.cohenD   = A.nodeStrength.meanDiff ./ (std(strDiff,0,1,'omitnan')+eps);
pNode = NaN(1, nG);
for k = 1:nG
    pNode(k) = pairedPvalue(strDiff(:,k), cfg.statsTest);
end
A.nodeStrength.pVals     = pNode;
A.nodeStrength.pAdj      = fdrBHlocal(pNode);
A.nodeStrength.signMask  = A.nodeStrength.pAdj < cfg.alpha;
end


%% ========================================================================
% =================  STATS PRIMITIVES  ===================================
%% ========================================================================

function p = pairedPvalue(d, testName)
d = d(~isnan(d));
if numel(d) < 3, p = NaN; return; end
switch lower(testName)
    case 'pairedt'
        m = mean(d); s = std(d); n = numel(d);
        if s == 0, p = double(m ~= 0); return; end
        t = m / (s/sqrt(n)); df = n - 1;
        p = betainc(df/(df + t^2), df/2, 0.5);
    case 'wilcoxon'
        d = d(d ~= 0); n = numel(d);
        if n < 5, p = NaN; return; end
        [~, ord] = sort(abs(d));
        ranks = zeros(1,n); ranks(ord) = 1:n;
        Wpos = sum(ranks(d > 0));
        mu = n*(n+1)/4; sg = sqrt(n*(n+1)*(2*n+1)/24);
        z = (Wpos - mu) / sg;
        p = erfc(abs(z)/sqrt(2));
    otherwise, error('Unknown statsTest: %s', testName);
end
end


function pAdj = fdrBHlocal(p)
p = p(:);
[ps, ix] = sort(p);
m = numel(p);
adj = ps .* m ./ (1:m)';
adj = min(1, cummin(adj, 'reverse'));
pAdj = nan(size(p));
pAdj(ix) = adj;
end


%% ========================================================================
% =================  CSV EXPORT  =========================================
%% ========================================================================

function writeScalarCSV(agg, cfg)
fid = fopen(fullfile(cfg.outputDir, 'global64_summary_stats.csv'), 'w');
fprintf(fid, ['metric,electrode_idx,electrode_label,mean_med,mean_non,', ...
              'mean_diff,sem_diff,cohens_d,p_value,p_adj_FDR,significant\n']);
metrics = {};
if cfg.runLZC, metrics{end+1} = 'LZC'; end
if cfg.runPAC, metrics{end+1} = 'PAC'; end
for mi = 1:numel(metrics)
    m = metrics{mi};
    A = agg.(m);
    for k = 1:numel(A.electrodes)
        e = A.electrodes(k);
        fprintf(fid, '%s,%d,%s,%.6f,%.6f,%.6f,%.6f,%.4f,%.6f,%.6f,%d\n', ...
                m, e, agg.refLabels{e}, ...
                A.meanMed(k), A.meanNon(k), A.meanDiff(k), ...
                A.semDiff(k), A.cohenD(k), A.pVals(k), A.pAdj(k), ...
                double(A.signMask(k)));
    end
end
fclose(fid);
fprintf('Wrote: global64_summary_stats.csv\n');
end


function writeMatrixCSV(agg, cfg)
A = agg.FC;
gE = A.electrodes;
nG = numel(gE);
fid = fopen(fullfile(cfg.outputDir, 'global64_FC_pairs.csv'), 'w');
fprintf(fid, ['elec_i_idx,elec_i_label,elec_j_idx,elec_j_label,', ...
              'mean_med,mean_non,mean_diff,cohens_d,p_value,p_adj_FDR,significant\n']);
for i = 1:nG-1
    for j = i+1:nG
        ei = gE(i); ej = gE(j);
        fprintf(fid, '%d,%s,%d,%s,%.6f,%.6f,%.6f,%.4f,%.6f,%.6f,%d\n', ...
                ei, agg.refLabels{ei}, ej, agg.refLabels{ej}, ...
                A.meanMed(i,j), A.meanNon(i,j), A.meanDiff(i,j), ...
                A.cohenD(i,j), A.pVals(i,j), A.pAdj(i,j), ...
                double(A.signMask(i,j)));
    end
end
fclose(fid);
fprintf('Wrote: global64_FC_pairs.csv\n');
end


%% ========================================================================
% =================  PLOTS - SCALAR METRICS  =============================
%% ========================================================================

function plotMetricTopoplots(agg, tag, prettyName, units, cfg, u)
A = agg.(tag);
gE = A.electrodes;

fullMed  = nan(numel(agg.refLabels), 1);
fullNon  = nan(numel(agg.refLabels), 1);
fullDiff = nan(numel(agg.refLabels), 1);
fullMed(gE)  = A.meanMed(:);
fullNon(gE)  = A.meanNon(:);
fullDiff(gE) = A.meanDiff(:);

fig = figure('Color','w','Position',[100 100 1400 460]);
clim_grp = safeClim([A.meanMed(:); A.meanNon(:)]);

subplot(1,3,1);
u.topoplotEEG(fullMed, agg.refXY, 'highlight', gE, 'cmap','parula', ...
              'clim', clim_grp, ...
              'title', sprintf('Meditators (n=%d pairs)', agg.nPairs));
cb = colorbar; cb.Label.String = units; cb.FontSize = 10;

subplot(1,3,2);
u.topoplotEEG(fullNon, agg.refXY, 'highlight', gE, 'cmap','parula', ...
              'clim', clim_grp, ...
              'title', sprintf('Non-Meditators (n=%d pairs)', agg.nPairs));
cb = colorbar; cb.Label.String = units; cb.FontSize = 10;

subplot(1,3,3);
clim_diff = safeAbsClim(A.meanDiff);
sigE = gE(A.signMask);
u.topoplotEEG(fullDiff, agg.refXY, 'highlight', sigE, ...
              'cmap', divergingCmap(), 'clim', clim_diff, ...
              'title', sprintf('Paired Delta  (%d/%d sig.)', numel(sigE), numel(gE)));
cb = colorbar; cb.Label.String = ['\Delta ' units]; cb.FontSize = 10;

sgtitle({sprintf('%s - 64-channel group comparison', prettyName), ...
         sprintf('Window: %s   |   FDR-corrected p < %.2f', ...
                 agg.windowLabel, cfg.alpha)}, ...
        'FontWeight', 'bold', 'FontSize', 13);

drawnow;
saveFigSafe(fig, fullfile(cfg.outputDir, sprintf('Fig_%s_topoplots', tag)));
end


function plotMetricSummary(agg, tag, prettyName, units, cfg, u)
A = agg.(tag);
gE = A.electrodes;
nG = numel(gE);
labels = agg.refLabels(gE);
cMed = [0.20 0.45 0.85];
cNon = [0.85 0.45 0.20];

fig = figure('Color','w','Position',[60 40 1500 1000]);

% Panel A: pooled-distribution KDE comparison.
subplot(2,2,1); hold on; box on; grid on;
medPool = A.medAll(:);
nonPool = A.nonAll(:);
[fM, xM] = ksdensitySimple(medPool);
[fN, xN] = ksdensitySimple(nonPool);
fill([xM xM(end:-1:1)], [fM zeros(1,numel(fM))], cMed, ...
     'EdgeColor', cMed*0.7, 'FaceAlpha', 0.45, 'LineWidth', 1.5, ...
     'DisplayName', 'Meditators');
fill([xN xN(end:-1:1)], [fN zeros(1,numel(fN))], cNon, ...
     'EdgeColor', cNon*0.7, 'FaceAlpha', 0.45, 'LineWidth', 1.5, ...
     'DisplayName', 'Non-Meditators');
mMed = mean(medPool,'omitnan');  mNon = mean(nonPool,'omitnan');
yL = ylim;
plot([mMed mMed], yL, '-', 'Color', cMed*0.5, 'LineWidth', 1.5, 'HandleVisibility','off');
plot([mNon mNon], yL, '-', 'Color', cNon*0.5, 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel(units); ylabel('Density');
title({'Pooled distribution across all 64 channels', ...
       sprintf('Med mean = %.4f   Non mean = %.4f', mMed, mNon)});
legend('Location','northeastoutside','Box','off');

% Panel B: scatter of per-electrode means.
subplot(2,2,2); hold on; box on; grid on;
mn = min([A.meanMed A.meanNon]);  mx = max([A.meanMed A.meanNon]);
plot([mn mx], [mn mx], '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
absD = abs(A.cohenD);
sz = 30 + 80 * (absD ./ max(absD + eps));
colors = mapToDiverging(A.cohenD);
for k = 1:nG
    plot(A.meanNon(k), A.meanMed(k), 'o', ...
         'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', ...
         'MarkerSize', sqrt(sz(k)));
    if A.signMask(k)
        text(A.meanNon(k), A.meanMed(k), [' ' labels{k}], 'FontSize', 8);
    end
end
xlabel(sprintf('Non-Meditator   %s', units));
ylabel(sprintf('Meditator   %s', units));
title({'Per-electrode means (paired)', ...
       'point color = Cohen''s d   |   labelled = FDR-sig.'});
axis equal; axis tight;

% Panel C: Cohen's d topoplot.
subplot(2,2,3);
fullD = nan(numel(agg.refLabels),1);
fullD(gE) = A.cohenD(:);
clim_d = safeAbsClim(A.cohenD);
u.topoplotEEG(fullD, agg.refXY, 'highlight', gE(A.signMask), ...
              'cmap', divergingCmap(), 'clim', clim_d, ...
              'title', 'Cohen''s d (paired effect size)');
cb = colorbar; cb.Label.String = 'Cohen''s d'; cb.FontSize = 10;

% Panel D: top-K differential electrodes paired bars.
subplot(2,2,4); hold on; box on; grid on;
K = min(cfg.topKElectrodes, nG);
[~, ord] = sort(abs(A.cohenD), 'descend');
ord = ord(1:K);
medK = A.medAll(:, ord);    nonK = A.nonAll(:, ord);
x = 1:K;  w = 0.30;
yM = mean(medK,1,'omitnan');  yN = mean(nonK,1,'omitnan');
bar(x - w/2, yM, w, 'FaceColor', cMed, 'EdgeColor','none', 'DisplayName', 'Meditator');
bar(x + w/2, yN, w, 'FaceColor', cNon, 'EdgeColor','none', 'DisplayName', 'Non-Meditator');
nP = size(medK, 1);
for p = 1:nP
    for k = 1:K
        if isnan(medK(p,k)) || isnan(nonK(p,k)), continue; end
        plot([x(k)-w/2, x(k)+w/2], [medK(p,k), nonK(p,k)], '-', ...
             'Color', [0.4 0.4 0.4 0.30], 'LineWidth', 0.6, ...
             'HandleVisibility', 'off');
    end
end
yL = ylim;
yStar = yL(2) + 0.03*diff(yL);
sigOrd = find(A.signMask(ord));
for kk = sigOrd(:)'
    text(x(kk), yStar, '*', 'HorizontalAlignment', 'center', ...
         'FontSize', 18, 'FontWeight', 'bold');
end
for k = 1:K
    text(x(k), yL(1)-0.07*diff(yL), sprintf('d=%+0.2f', A.cohenD(ord(k))), ...
         'HorizontalAlignment','center','FontSize',8,'Color',[0.3 0.3 0.3]);
end
ylim([yL(1)-0.12*diff(yL), yL(2)+0.10*diff(yL)]);
set(gca, 'XTick', x, 'XTickLabel', labels(ord), 'XTickLabelRotation', 45);
ylabel(units);
title(sprintf('Top %d electrodes by |Cohen''s d|  (* = FDR p<%.2f)', K, cfg.alpha));
legend('Location','northeastoutside','Box','off');

sgtitle({sprintf('%s - 64-channel summary', prettyName), ...
         sprintf('Pairs = %d   |   sig elecs = %d', agg.nPairs, sum(A.signMask))}, ...
        'FontWeight', 'bold', 'FontSize', 13);

drawnow;
saveFigSafe(fig, fullfile(cfg.outputDir, sprintf('Fig_%s_summary', tag)));
end


%% ========================================================================
% =================  PLOTS - FUNCTIONAL CONNECTIVITY  ====================
%% ========================================================================

function plotFCMatrices(agg, cfg, u)
A = agg.FC;
gE = A.electrodes;
labels = agg.refLabels(gE);
xy = agg.refXY(gE, :);
[~, ord] = sortrows([-xy(:,2), xy(:,1)]);
labOrd = labels(ord);
nG = numel(ord);

fig = figure('Color','w','Position',[40 40 1600 540]);

% Use percentile-based clim so a few extreme values don't wash out the rest.
allVals = [A.meanMed(:); A.meanNon(:)];
allVals = allVals(~isnan(allVals) & isfinite(allVals));
if isempty(allVals)
    clim_grp = [0 1];
else
    clim_grp = [0, prctile(allVals, cfg.fcMatrixPctile)];
    if clim_grp(2) <= 0, clim_grp(2) = max(allVals); end
    if clim_grp(2) <= clim_grp(1), clim_grp = [0 max(eps, max(allVals))]; end
end

% Tick spacing - with 64 electrodes, label every 4th electrode for clarity.
tickStride = max(1, round(nG/16));
tickPos = 1:tickStride:nG;
tickLab = labOrd(tickPos);

subplot(1,3,1);
imagesc(A.meanMed(ord,ord), clim_grp); axis square;
colormap(gca,'parula');
applyMatrixAxes(gca, tickPos, tickLab);
cb = colorbar; cb.Label.String = 'Mutual Information (bits)'; cb.FontSize = 10;
title(sprintf('Meditators (n=%d pairs)', agg.nPairs), 'FontSize', 12);

subplot(1,3,2);
imagesc(A.meanNon(ord,ord), clim_grp); axis square;
colormap(gca,'parula');
applyMatrixAxes(gca, tickPos, tickLab);
cb = colorbar; cb.Label.String = 'Mutual Information (bits)'; cb.FontSize = 10;
title(sprintf('Non-Meditators (n=%d pairs)', agg.nPairs), 'FontSize', 12);

subplot(1,3,3);
D = A.meanDiff(ord,ord);
clim_diff = safeAbsClim(D);
imagesc(D, clim_diff); axis square;
colormap(gca, divergingCmap()); cb = colorbar; cb.Label.String = '\Delta MI (bits)'; cb.FontSize = 10;
hold on;
sigOrd = A.signMask(ord, ord);
[ii, jj] = find(triu(sigOrd, 1));
plot(jj, ii, 'k.', 'MarkerSize', 6, 'HandleVisibility','off');
plot(ii, jj, 'k.', 'MarkerSize', 6, 'HandleVisibility','off');
applyMatrixAxes(gca, tickPos, tickLab);
title({'Paired Delta  (Med - Non-Med)', sprintf('dots = FDR p<%.2f', cfg.alpha)}, ...
      'FontSize', 12);

sgtitle({'Functional Connectivity (Mutual Information) - 64-channel comparison', ...
         sprintf('Window: %s   |   electrodes ordered anterior-to-posterior', ...
                 agg.windowLabel)}, ...
        'FontWeight', 'bold', 'FontSize', 13);

drawnow;
saveFigSafe(fig, fullfile(cfg.outputDir, 'Fig_FC_matrices'));
end


function applyMatrixAxes(ax, tickPos, tickLab)
set(ax, 'XTick', tickPos, 'XTickLabel', tickLab, ...
        'YTick', tickPos, 'YTickLabel', tickLab, ...
        'XTickLabelRotation', 60, 'FontSize', 8, ...
        'TickLength', [0.005 0.005]);
xlabel(''); ylabel('');
end


function plotFCSummary(agg, cfg, u)
A = agg.FC;
gE = A.electrodes;
labels = agg.refLabels(gE);
nG = numel(gE);

[ii, jj] = find(triu(true(nG), 1));
nE = numel(ii);
diffVec  = arrayfun(@(k) A.meanDiff(ii(k), jj(k)), 1:nE);
dVec     = arrayfun(@(k) A.cohenD  (ii(k), jj(k)), 1:nE);
sigVec   = arrayfun(@(k) A.signMask(ii(k), jj(k)), 1:nE);

[~, ord] = sort(abs(dVec), 'descend');
K = min(cfg.topKConnections, nE);
topIdx = ord(1:K);

fig = figure('Color','w','Position',[40 40 1500 1000]);

% Panel A: top-K connections horizontal bar.
subplot(2,2,1); hold on; box on;
yvals = (1:K)';
for k = 1:K
    idx = topIdx(k);
    if diffVec(idx) >= 0, col = [0.85 0.20 0.20]; else, col = [0.20 0.40 0.85]; end
    barh(yvals(k), diffVec(idx), 'FaceColor', col, 'EdgeColor', 'none');
    if sigVec(idx)
        text(diffVec(idx), yvals(k), '*', 'FontSize', 14, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'left');
    end
end
ytl = arrayfun(@(k) sprintf('%s - %s  (d=%+0.2f)', ...
       labels{ii(topIdx(k))}, labels{jj(topIdx(k))}, dVec(topIdx(k))), ...
       (1:K)', 'UniformOutput', false);
set(gca,'YTick', yvals, 'YTickLabel', ytl, 'YDir','reverse', 'FontSize', 8);
xline(0,'k-');
xlabel('\Delta Mutual Information (bits)   [Med - Non-Med]');
title(sprintf('Top %d differential connections  (* = FDR p<%.2f)', K, cfg.alpha));

% Panel B: most-affected nodes.
subplot(2,2,2); hold on; box on; grid on;
nodeImpact = sum(abs(A.cohenD), 2, 'omitnan');
[nodeImpact_s, nord] = sort(nodeImpact, 'descend');
K2 = min(nG, 15);
barh(1:K2, nodeImpact_s(1:K2), 'FaceColor', [0.45 0.55 0.70], 'EdgeColor','none');
set(gca,'YTick', 1:K2, 'YTickLabel', labels(nord(1:K2)), 'YDir','reverse', 'FontSize', 9);
xlabel('Sum of |Cohen''s d| across all connections');
title({'Most differentially-connected channels', ...
       'electrodes whose connections differ most across groups'});

% Panel C: node-strength Delta topoplot.
subplot(2,2,3);
NS = A.nodeStrength;
fullD = nan(numel(agg.refLabels),1);
fullD(gE) = NS.meanDiff(:);
clim_d = safeAbsClim(NS.meanDiff);
u.topoplotEEG(fullD, agg.refXY, 'highlight', gE(NS.signMask), ...
              'cmap', divergingCmap(), 'clim', clim_d, ...
              'title', 'Node strength Delta  (Med - Non-Med)');
cb = colorbar; cb.Label.String = '\Delta\Sigma MI'; cb.FontSize = 10;

% Panel D: per-connection scatter.
subplot(2,2,4); hold on; box on; grid on;
allMed = arrayfun(@(k) A.meanMed(ii(k), jj(k)), 1:nE);
allNon = arrayfun(@(k) A.meanNon(ii(k), jj(k)), 1:nE);
mn = min([allMed allNon]); mx = max([allMed allNon]);
plot([mn mx], [mn mx], '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
nonsig = ~sigVec;
plot(allNon(nonsig), allMed(nonsig), 'o', ...
     'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerSize', 4, ...
     'DisplayName', sprintf('non-sig (n=%d)', sum(nonsig)));
sigPos = sigVec & diffVec > 0;
sigNeg = sigVec & diffVec < 0;
plot(allNon(sigPos), allMed(sigPos), 'o', ...
     'MarkerFaceColor', [0.85 0.20 0.20], 'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
     'DisplayName', sprintf('Med > Non (n=%d)', sum(sigPos)));
plot(allNon(sigNeg), allMed(sigNeg), 'o', ...
     'MarkerFaceColor', [0.20 0.40 0.85], 'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
     'DisplayName', sprintf('Non > Med (n=%d)', sum(sigNeg)));
xlabel('Non-Meditator MI (bits)'); ylabel('Meditator MI (bits)');
axis equal; axis tight;
legend('Location','northeastoutside','Box','off');
title('Per-connection means  (FDR-significant highlighted)');

sgtitle({'Functional Connectivity (Mutual Information) - summary', ...
         sprintf('Pairs = %d   |   sig connections = %d / %d', ...
                 agg.nPairs, sum(sigVec), nE)}, ...
        'FontWeight', 'bold', 'FontSize', 13);

drawnow;
saveFigSafe(fig, fullfile(cfg.outputDir, 'Fig_FC_summary'));
end


function plotFCNetwork(agg, cfg, u)
A = agg.FC;
gE = A.electrodes;
xy = agg.refXY;

fig = figure('Color','w','Position',[40 40 1500 520]);

subplot(1,3,1);
drawHeadOutline();
drawTopEdgesLocal(A.meanMed, gE, xy, 0.10, [0.20 0.45 0.85]);
drawNodesLocal(xy, agg.refLabels, gE);
title(sprintf('Meditators - top 10%% edges (n=%d pairs)', agg.nPairs), 'FontSize', 12);

subplot(1,3,2);
drawHeadOutline();
drawTopEdgesLocal(A.meanNon, gE, xy, 0.10, [0.85 0.45 0.20]);
drawNodesLocal(xy, agg.refLabels, gE);
title(sprintf('Non-Meditators - top 10%% edges (n=%d pairs)', agg.nPairs), 'FontSize', 12);

subplot(1,3,3);
drawHeadOutline();
nG = numel(gE);
Dsig = zeros(nG, nG);
Dsig(A.signMask) = A.meanDiff(A.signMask);
drawDiffEdgesLocal(Dsig, gE, xy);
drawNodesLocal(xy, agg.refLabels, gE);
title({'FDR-significant Delta  edges', 'red: Med>Non   blue: Non>Med'}, 'FontSize', 12);

sgtitle({'Functional Connectivity (Mutual Information) Network', ...
         sprintf('Window: %s   |   alpha=%.2f', agg.windowLabel, cfg.alpha)}, ...
        'FontWeight', 'bold', 'FontSize', 13);

drawnow;
saveFigSafe(fig, fullfile(cfg.outputDir, 'Fig_FC_network'));
end


%% ========================================================================
% =================  GLOBAL SUMMARY (one figure)  ========================
%% ========================================================================

function plotGlobalSummary(agg, cfg, u)
have = [cfg.runLZC, cfg.runPAC, cfg.runFC];
nMetrics = sum(have);
if nMetrics == 0, return; end

fig = figure('Color','w','Position',[40 40 480*nMetrics 820]);

col = 0;
if cfg.runLZC
    col = col + 1;
    plotSummaryColumn(agg.LZC, agg, 'LZC', 'normalised LZC', col, nMetrics, u);
end
if cfg.runPAC
    col = col + 1;
    plotSummaryColumn(agg.PAC, agg, 'PAC (Modulation Index)', ...
                      'Modulation Index', col, nMetrics, u);
end
if cfg.runFC
    col = col + 1;
    A = agg.FC.nodeStrength;
    A.electrodes = agg.FC.electrodes;
    plotSummaryColumn(A, agg, 'FC node strength (Mutual Information)', ...
                      '\Sigma MI (bits)', col, nMetrics, u);
end

sgtitle({'Global Summary - meditator vs non-meditator differences (all 64 channels)', ...
         sprintf('Pairs = %d   |   FDR p<%.2f', agg.nPairs, cfg.alpha)}, ...
        'FontWeight', 'bold', 'FontSize', 13);

drawnow;
saveFigSafe(fig, fullfile(cfg.outputDir, 'Fig_GlobalSummary'));
end


function plotSummaryColumn(A, agg, name, units, col, nCol, u)
gE = A.electrodes;

subplot(2, nCol, col);
fullD = nan(numel(agg.refLabels), 1);
fullD(gE) = A.meanDiff(:);
clim_diff = safeAbsClim(A.meanDiff);
sigE = gE(A.signMask);
u.topoplotEEG(fullD, agg.refXY, 'highlight', sigE, ...
              'cmap', divergingCmap(), 'clim', clim_diff, ...
              'title', sprintf('%s - Delta  (Med - Non)', name));
cb = colorbar; cb.Label.String = ['\Delta ' units]; cb.FontSize = 10;

subplot(2, nCol, col + nCol);
fullCd = nan(numel(agg.refLabels), 1);
fullCd(gE) = A.cohenD(:);
clim_d = safeAbsClim(A.cohenD);
u.topoplotEEG(fullCd, agg.refXY, 'highlight', sigE, ...
              'cmap', divergingCmap(), 'clim', clim_d, ...
              'title', sprintf('%s - Cohen''s d', name));
cb = colorbar; cb.Label.String = 'Cohen''s d'; cb.FontSize = 10;
end


%% ========================================================================
% =================  PER-PAIR FIGURE  ====================================
%% ========================================================================

function plotPairFigure(R, agg, p, nP, pairsDir, cfg, u)
gE = agg.globalElecs;

fig = figure('Color','w','Position',[40 40 1500 1300]);

if cfg.runLZC
    fullM = nan(numel(agg.refLabels),1);  fullM(gE) = R.medLZC(gE);
    fullN = nan(numel(agg.refLabels),1);  fullN(gE) = R.nonLZC(gE);
    diffV = R.medLZC(gE) - R.nonLZC(gE);
    fullD = nan(numel(agg.refLabels),1);  fullD(gE) = diffV;
    cl     = safeClim([fullM; fullN]);
    clDiff = safeAbsClim(diffV);
    subplot(3,3,1); u.topoplotEEG(fullM, agg.refXY, 'highlight', gE, ...
        'cmap','parula','clim',cl, 'title', sprintf('LZC - Med (%s)', R.medID));
    cb = colorbar; cb.Label.String = 'LZC';
    subplot(3,3,2); u.topoplotEEG(fullN, agg.refXY, 'highlight', gE, ...
        'cmap','parula','clim',cl, 'title', sprintf('LZC - Non (%s)', R.nonID));
    cb = colorbar; cb.Label.String = 'LZC';
    subplot(3,3,3);
    u.topoplotEEG(fullD, agg.refXY, 'highlight', gE, ...
        'cmap', divergingCmap(),'clim',clDiff, 'title', 'LZC - Delta ');
    cb = colorbar; cb.Label.String = '\Delta LZC';
end

if cfg.runPAC
    fullM = nan(numel(agg.refLabels),1);  fullM(gE) = R.medPAC(gE);
    fullN = nan(numel(agg.refLabels),1);  fullN(gE) = R.nonPAC(gE);
    diffV = R.medPAC(gE) - R.nonPAC(gE);
    fullD = nan(numel(agg.refLabels),1);  fullD(gE) = diffV;
    cl     = safeClim([fullM; fullN]);
    clDiff = safeAbsClim(diffV);
    subplot(3,3,4); u.topoplotEEG(fullM, agg.refXY, 'highlight', gE, ...
        'cmap','parula','clim',cl, ...
        'title', sprintf('PAC (Modulation Index) - Med (%s)', R.medID));
    cb = colorbar; cb.Label.String = 'Modulation Index';
    subplot(3,3,5); u.topoplotEEG(fullN, agg.refXY, 'highlight', gE, ...
        'cmap','parula','clim',cl, ...
        'title', sprintf('PAC (Modulation Index) - Non (%s)', R.nonID));
    cb = colorbar; cb.Label.String = 'Modulation Index';
    subplot(3,3,6);
    u.topoplotEEG(fullD, agg.refXY, 'highlight', gE, ...
        'cmap', divergingCmap(),'clim',clDiff, ...
        'title', 'PAC (Modulation Index) - Delta ');
    cb = colorbar; cb.Label.String = '\Delta MI';
end

if cfg.runFC && numel(gE) >= 2
    Mm = R.medFC(gE, gE);  Mn = R.nonFC(gE, gE);  Md = Mm - Mn;
    xy = agg.refXY(gE,:);
    [~, ord] = sortrows([-xy(:,2), xy(:,1)]);

    % Percentile-based clim for the Med/Non panels.
    allVals = [Mm(:); Mn(:)];
    allVals = allVals(~isnan(allVals) & isfinite(allVals));
    if isempty(allVals)
        cl = [0 1];
    else
        cl = [0, prctile(allVals, cfg.fcMatrixPctile)];
        if cl(2) <= cl(1), cl = [0 max(eps, max(allVals))]; end
    end
    clDiff = safeAbsClim(Md);

    subplot(3,3,7);
    imagesc(Mm(ord,ord), cl); axis square; colormap(gca,'parula');
    set(gca,'XTick',[],'YTick',[]); cb = colorbar; cb.Label.String = 'Mutual Information (bits)';
    title(sprintf('FC (Mutual Information) - Med (%s)', R.medID));

    subplot(3,3,8);
    imagesc(Mn(ord,ord), cl); axis square; colormap(gca,'parula');
    set(gca,'XTick',[],'YTick',[]); cb = colorbar; cb.Label.String = 'Mutual Information (bits)';
    title(sprintf('FC (Mutual Information) - Non (%s)', R.nonID));

    subplot(3,3,9);
    imagesc(Md(ord,ord), clDiff); axis square; colormap(gca, divergingCmap());
    set(gca,'XTick',[],'YTick',[]); cb = colorbar; cb.Label.String = '\Delta MI (bits)';
    title('FC (Mutual Information) - Delta ');
end

sgtitle({sprintf('Pair %d/%d - %s vs %s', p, nP, R.medID, R.nonID), ...
         sprintf('Window: %s   |   64 channels', R.windowLabel)}, ...
        'FontWeight', 'bold', 'FontSize', 13);

drawnow;
baseFile = fullfile(pairsDir, sprintf('pair_%s_%s', R.medID, R.nonID));
saveFigSafe(fig, baseFile);
fprintf('   [%2d/%2d] %s vs %s\n', p, nP, R.medID, R.nonID);
end


%% ========================================================================
% =================  FIGURE SAVE WRAPPER  ================================
%% ========================================================================

function saveFigSafe(fig, baseFile)
% SAVEFIGSAFE  Save figure as PNG (300 dpi) and editable .fig.
%
%   The figure remains visible on screen after saving so that it can be
%   inspected or manually exported elsewhere.  drawnow() is called first
%   to flush all pending graphics state before serialisation - this
%   prevents the "loadable .png but unloadable .fig" failure mode that
%   occurs when MATLAB's graphics state lags behind the figure object
%   tree.
set(fig, 'PaperPositionMode', 'auto');
drawnow;
print(fig, [baseFile '.png'], '-dpng', '-r300');
savefig(fig, [baseFile '.fig'], 'compact');
end


%% ========================================================================
% =================  PLOT PRIMITIVES  ====================================
%% ========================================================================

function colors = mapToDiverging(vals)
v = vals(:);
mx = max(abs(v), [], 'omitnan');
if mx <= 0 || isnan(mx), mx = eps; end
v = v / mx;
cmap = divergingCmap();
n = size(cmap, 1);
idx = round((v + 1) / 2 * (n - 1)) + 1;
idx(isnan(idx)) = round((n+1)/2);
idx = max(1, min(n, idx));
colors = cmap(idx, :);
end


function cmap = divergingCmap()
% RdBu-style diverging colormap, blue (negative) -> white -> red (positive).
n = 128;
b = [linspace(0.10,1,n)' linspace(0.30,1,n)' linspace(0.85,1,n)'];
r = [linspace(1,0.85,n)' linspace(1,0.15,n)' linspace(1,0.15,n)'];
cmap = [b; r];
end


function cl = safeClim(vals, fallback)
if nargin < 2, fallback = [0 1]; end
v = vals(:);
v = v(~isnan(v) & isfinite(v));
if isempty(v) || min(v) >= max(v)
    cl = fallback;
else
    cl = [min(v) max(v)];
end
end


function cl = safeAbsClim(vals, fallback)
if nargin < 2, fallback = eps; end
v = vals(:);
v = v(~isnan(v) & isfinite(v));
if isempty(v)
    mx = fallback;
else
    mx = max(abs(v));
    if mx <= 0, mx = fallback; end
end
cl = [-mx, mx];
end


function drawHeadOutline()
hold on; axis equal; axis off;
xlim([-0.6 0.6]); ylim([-0.6 0.6]);
th = linspace(0, 2*pi, 200);
plot(0.5*cos(th), 0.5*sin(th), 'k', 'LineWidth', 1.4);
plot([-0.05 0 0.05], [0.5 0.55 0.5], 'k', 'LineWidth', 1.4);
plot([-0.5 -0.55 -0.55 -0.5], [0.1 0.05 -0.05 -0.1], 'k', 'LineWidth', 1.2);
plot([ 0.5  0.55  0.55  0.5], [0.1 0.05 -0.05 -0.1], 'k', 'LineWidth', 1.2);
end


function drawTopEdgesLocal(M, elecIdx, xy, frac, edgeColor)
% Draw the top fraction of edges among the selected channels.
%   M       : [nG x nG] connectivity matrix in LOCAL global-electrode basis
%   elecIdx : [1 x nG] mapping from local index to original electrode index
%   xy      : [Neeg x 2] electrode positions in original electrode basis
nG = numel(elecIdx);
[ii, jj] = find(triu(true(nG), 1));
vals = arrayfun(@(p) M(ii(p), jj(p)), 1:numel(ii));
ok = ~isnan(vals);
ii = ii(ok); jj = jj(ok); vals = vals(ok);
if isempty(vals), return; end
[~, sIdx] = sort(vals, 'descend');
nKeep = max(1, round(frac * numel(vals)));
sIdx = sIdx(1:nKeep);
vMin = min(vals(sIdx)); vMax = max(vals(sIdx));
if vMax == vMin, vMax = vMin + eps; end
for k = 1:numel(sIdx)
    pp = sIdx(k);
    a = elecIdx(ii(pp)); b = elecIdx(jj(pp));
    if any(isnan(xy(a,:))) || any(isnan(xy(b,:))), continue; end
    w   = (vals(pp) - vMin) / (vMax - vMin);
    lw  = 0.4 + 2.4 * w;
    alp = 0.25 + 0.65 * w;
    plot(xy([a b],1), xy([a b],2), '-', ...
         'Color', [edgeColor alp], 'LineWidth', lw);
end
end


function drawDiffEdgesLocal(Dsig, elecIdx, xy)
nG = numel(elecIdx);
[ii, jj] = find(triu(true(nG), 1));
vals = arrayfun(@(p) Dsig(ii(p), jj(p)), 1:numel(ii));
ok = vals ~= 0 & ~isnan(vals);
ii = ii(ok); jj = jj(ok); vals = vals(ok);
if isempty(vals), return; end
aMax = max(abs(vals));
if aMax == 0, aMax = eps; end
for k = 1:numel(vals)
    a = elecIdx(ii(k)); b = elecIdx(jj(k));
    if any(isnan(xy(a,:))) || any(isnan(xy(b,:))), continue; end
    w  = abs(vals(k)) / aMax;
    lw = 0.6 + 2.6 * w;
    al = 0.40 + 0.55 * w;
    if vals(k) > 0
        col = [0.85 0.15 0.15 al];
    else
        col = [0.15 0.30 0.85 al];
    end
    plot(xy([a b],1), xy([a b],2), '-', 'Color', col, 'LineWidth', lw);
end
end


function drawNodesLocal(xy, eegLabels, elecIdx)
keep = ~any(isnan(xy), 2);
plot(xy(keep,1), xy(keep,2), 'k.', 'MarkerSize', 7);
foc = elecIdx(:);
foc = foc(~any(isnan(xy(foc,:)), 2));
if numel(foc) > 30
    plot(xy(foc,1), xy(foc,2), 'ko', ...
         'MarkerFaceColor','w','MarkerSize',5,'LineWidth',0.8);
else
    if ~isempty(foc)
        plot(xy(foc,1), xy(foc,2), 'ko', ...
             'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.2);
        for k = 1:numel(foc)
            text(xy(foc(k),1), xy(foc(k),2)+0.035, eegLabels{foc(k)}, ...
                 'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
        end
    end
end
end


function [f, x] = ksdensitySimple(v)
v = v(~isnan(v));
if isempty(v), f = 0; x = 0; return; end
n = numel(v);
sd = std(v); if sd == 0, sd = 1; end
h = 1.06 * sd * n^(-1/5);
if h <= 0, h = 0.1; end
xmin = min(v) - 3*h;  xmax = max(v) + 3*h;
x = linspace(xmin, xmax, 100);
f = zeros(1, numel(x));
for i = 1:numel(x)
    f(i) = sum(exp(-((x(i) - v).^2) / (2*h^2))) / (n * h * sqrt(2*pi));
end
end


%% ========================================================================
% =================  MATH PRIMITIVES  ====================================
%% ========================================================================

function c = lz76local(s)
if ~ischar(s), s = char(s(:)' + '0'); end
n = length(s);
if n == 0, c = 0; return; end
c = 1; S = s(1); Q = '';
for i = 2:n
    Q = [Q s(i)]; %#ok<AGROW>
    SQv = [S Q(1:end-1)];
    if isempty(strfind(SQv, Q)) %#ok<STREMP>
        S = [S Q]; %#ok<AGROW>
        Q = '';
        c = c + 1;
    end
end
if ~isempty(Q), c = c + 1; end
end


function [phs, amp] = trialwiseHilbertLocal(phaseSig, ampSig, tIdx)
nT = size(phaseSig,1);
phsList = cell(nT,1); ampList = cell(nT,1);
for ti = 1:nT
    p = angle(hilbert(phaseSig(ti,:)));
    a = abs  (hilbert(ampSig  (ti,:)));
    phsList{ti} = p(tIdx);
    ampList{ti} = a(tIdx);
end
phs = cat(2, phsList{:})';
amp = cat(2, ampList{:})';
end


function MI = tortMIlocal(phs, amp, nbins)
% Tort 2010 modulation index: phase histogram of mean amplitude.
edges = linspace(-pi, pi, nbins+1);
[~, ~, bin] = histcounts(phs, edges);
meanAmp = zeros(nbins,1);
for k = 1:nbins
    sel = (bin == k);
    if any(sel), meanAmp(k) = mean(amp(sel)); end
end
if sum(meanAmp) <= 0, MI = NaN; return; end
P = meanAmp / sum(meanAmp);
H = -sum(P .* log(P + eps));
MI = (log(nbins) - H) / log(nbins);
end


function y = bandpassRowsLocal(x, Fs, lo, hi)
order = round(Fs/2);
if hi >= Fs/2, hi = Fs/2 - 1; end
if lo <= 0,    lo = 0.5;       end
b = fir1(order, [lo hi]/(Fs/2), 'bandpass');
y = zeros(size(x));
for r = 1:size(x,1)
    y(r,:) = filtfilt(b, 1, x(r,:));
end
end


function bx = quantileBinLocal(x, nbins)
edges = quantile(x(:), linspace(0,1,nbins+1));
edges = unique(edges);
if numel(edges) < 2
    bx = ones(1, numel(x), 'uint8'); return;
end
edges(1)   = -inf;
edges(end) =  inf;
[~, ~, bin] = histcounts(x, edges);
bin(bin == 0) = 1;
bx = uint8(bin(:)');
end


function mi = miFromBinsLocal(bx, by, nbins)
bx = double(bx(:));  by = double(by(:));
J = accumarray([bx by], 1, [nbins nbins]);
N = sum(J(:));
if N == 0, mi = NaN; return; end
J = J / N;
px = sum(J,2);  py = sum(J,1);
PxPy = px * py;
mask = (J > 0) & (PxPy > 0);
mi = sum(J(mask) .* log2(J(mask) ./ PxPy(mask)));
end


%% ========================================================================
% =================  PATH / IO HELPERS  ==================================
%% ========================================================================

function lfp = findLFPPath(rootDir, subjID, protocol, dateOverride)
lfp = '';
candidates = {};
candidates{end+1} = fullfile(rootDir, subjID, protocol);
candidates{end+1} = fullfile(rootDir, subjID, protocol, 'LFP');
candidates{end+1} = fullfile(rootDir, subjID, protocol, 'segmentedData', 'LFP');
eegDir = fullfile(rootDir, subjID, 'EEG');
if exist(eegDir, 'dir')
    if nargin >= 4 && ~isempty(dateOverride)
        candidates{end+1} = fullfile(eegDir, dateOverride, protocol, ...
                                     'segmentedData', 'LFP');
    else
        d = dir(eegDir);
        d = d([d.isdir] & ~startsWith({d.name}, '.'));
        if ~isempty(d)
            [~, ord] = sort({d.name});
            d = d(ord(end:-1:1));
            for k = 1:numel(d)
                candidates{end+1} = fullfile(eegDir, d(k).name, protocol, ...
                                             'segmentedData', 'LFP'); %#ok<AGROW>
            end
        end
    end
end
for k = 1:numel(candidates)
    c = candidates{k};
    if exist(c, 'dir') && exist(fullfile(c, 'badTrials_wo_v8.mat'), 'file')
        lfp = c; return;
    end
end
end


function pairs = readPairList(csvPath)
pairs = {};
fid = fopen(csvPath, 'r');
if fid < 0, return; end
sawHeader = false;
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || startsWith(line, '#'), continue; end
    cols = strsplit(line, ',');
    cols = cellfun(@strtrim, cols, 'UniformOutput', false);
    while numel(cols) < 4, cols{end+1} = ''; end %#ok<AGROW>
    cols = cols(1:4);
    if ~sawHeader && any(~cellfun(@isempty, regexpi(cols, ...
            '^(meditator|control|med|non|subject)', 'match', 'once')))
        sawHeader = true; continue;
    end
    sawHeader = true;
    pairs(end+1, 1:4) = cols; %#ok<AGROW>
end
fclose(fid);
end


function [medID, nonID, medDate, nonDate] = unpackPairRow(pairs, p)
medID = pairs{p,1}; nonID = pairs{p,2};
medDate = ''; nonDate = '';
if size(pairs,2) >= 3, medDate = pairs{p,3}; end
if size(pairs,2) >= 4, nonDate = pairs{p,4}; end
end


function elecs = pickGoodElecs(S, cfg)
% Decide which electrodes the per-subject compute helpers operate on.
%   useAllElectrodes = true  -> every channel that has a data file
%                               (per-electrode QC is bypassed)
%   useAllElectrodes = false -> data-file-present AND not in
%                               the subject's flagged-bad list
if isfield(cfg, 'useAllElectrodes') && cfg.useAllElectrodes
    elecs = find(S.elecAvailable);
else
    elecs = setdiff(find(S.elecAvailable), S.badElectrodes);
end
end


function mkdirSafe(d)
if ~exist(d, 'dir'), mkdir(d); end
end
