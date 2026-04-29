function utils = eegUtils()
% EEGUTILS  Shared helpers for the Meditators vs Non-Meditators EEG project.
% =========================================================================
%
%   PURPOSE
%   -------
%   Shared helpers used by all three analysis scripts (Lempel-Ziv,
%   Phase-Amplitude Coupling, Mutual Information).  Keeping common
%   routines here so that every script will:
%       (a) load subject data with identical conventions,
%       (b) match each EEG sample to the correct moment in time using the
%           timeVals vector that was saved with the recording,
%       (c) reject the same bad trials & bad electrodes,
%       (d) lay out topoplots with the same coordinate convention,
%       (e) save figures the same way (high-res .png + editable .fig).
%
%   USAGE
%   -----
%       u = eegUtils();
%       S = u.loadSubject(folder);          % load one subject's data
%       [tIdx, t] = u.getTimeWindowIdx(S.timeVals, [-1 1]);   % window
%       y = u.bandpassFilter(x, S.Fs, 8, 13);                 % filter
%       u.topoplotEEG(values, S.elecXY, 'title','Alpha power');
%
%   The function returns a struct of function handles.  Treat 'u' as a
%   namespace.
%
%   Author: Kirubananth S, EEG Meditator project, Neural Signal Processing Course 2026.
% =========================================================================

% Public handles.  Add new helpers below in their own local function.
utils.loadSubject            = @loadSubject;
utils.getTimeWindowIdx       = @getTimeWindowIdx;
utils.getEEGPositions        = @getEEGPositions;
utils.bandpassFilter         = @bandpassFilter;
utils.topoplotEEG            = @topoplotEEG;
utils.saveFigDual            = @saveFigDual;
utils.fdrBH                  = @fdrBH;
utils.defaultFrontoParietal  = @defaultFrontoParietal;
utils.printDivider           = @printDivider;
utils.formatTimeWindow       = @formatTimeWindow;
end


%% ========================================================================
function S = loadSubject(folder)
% LOADSUBJECT  Read all metadata + per-electrode data for one subject.
% -------------------------------------------------------------------------
%   The subject folder is expected to contain:
%       lfpInfo.mat            - timeVals, Fs (sampling rate), electrode list
%       badTrials_wo_v8.mat    - badTrials, badElecs, electrode labels,
%                                highPriorityElectrodeList
%       impedanceData.mat      - (optional) per-electrode impedance
%       elec1.mat ... elecN.mat - one .mat per electrode containing
%                                 analogData [nTrials x nSamples]
%                                 analogInfo (theta, radius, X, Y, Z)
%
%   The returned struct S has these fields used downstream:
%       S.folder                 (char)        folder path passed in
%       S.timeVals               (1 x nSamp)   time-axis in seconds, the
%                                              SAME vector used to record
%                                              every elec*.mat file
%       S.Fs                     (scalar)      sampling rate (Hz),
%                                              derived from timeVals
%       S.totalTrials            (scalar)      number of trials per elec
%       S.badTrials              (col vec)     trial indices to drop
%       S.eegLabels              (Neeg x 1)    cell array of '10-10' labels
%       S.highPriorityElecs      (1 x ?)       indices, occipito-parietal
%       S.badElectrodes          (1 x ?)       electrodes to drop (union of
%                                              noisyElecs + flatPSDElecs +
%                                              badImpedanceElecs +
%                                              declaredBadElectrodes)
%       S.electrodeImpedances    (col vec)     impedance per electrode
%       S.data                   (Neeg x nTr x nSamp, single)
%                                              EEG data; NaN for missing
%                                              electrodes
%       S.elecAvailable          (1 x Neeg)    logical, true if elec*.mat
%                                              was found for this elec
%       S.elecXY                 (Neeg x 2)    2-D topoplot coordinates
%                                              (preferred from analogInfo,
%                                              fallback to standard 10-10)
% -------------------------------------------------------------------------

    S = struct();
    S.folder = folder;

    % --------------------------------------------------------------- %
    %  1. lfpInfo  -- defines the time axis used by EVERY electrode.   %
    % --------------------------------------------------------------- %
    li = load(fullfile(folder,'lfpInfo.mat'));
    S.timeVals          = li.timeVals(:)';
    % Sampling rate inferred from the time vector itself.  Using
    % median(diff()) is robust against any spurious sample at the edges.
    S.Fs                = round(1/median(diff(S.timeVals)));
    S.electrodesStored  = double(li.electrodesStored(:)');
    S.analogChannelsStored = double(li.analogChannelsStored(:)');

    % A 1000 Hz file should have ~1 ms steps - warn if it doesn't.
    if abs(S.Fs - 1/median(diff(S.timeVals))) > 0.5
        warning('eegUtils:FsRoundingMismatch', ...
            'Fs rounded to %d Hz but median diff(timeVals) was %.6f s', ...
            S.Fs, median(diff(S.timeVals)));
    end

    % --------------------------------------------------------------- %
    %  2. bad trials, bad electrodes, channel labels                    %
    % --------------------------------------------------------------- %
    bt = load(fullfile(folder,'badTrials_wo_v8.mat'));
    S.badTrials          = double(bt.badTrials(:));
    S.totalTrials        = double(bt.totalTrials);
    S.highPriorityElecs  = double(bt.highPriorityElectrodeList(:)');

    % Convert the cellstr labels to a clean column cell of chars.
    rawLabels = bt.eegElectrodeLabels(:);
    S.eegLabels = cell(numel(rawLabels),1);
    for k = 1:numel(rawLabels)
        v = rawLabels{k};
        if iscell(v), v = v{1}; end
        S.eegLabels{k} = char(v);
    end
    Neeg = numel(S.eegLabels);

    % Union of every recognized "bad" sub-list inside the badElecs struct.
    % These electrodes are dropped from analyses subject-by-subject.
    badE = [];
    if isfield(bt,'badElecs') && ~isempty(bt.badElecs)
        be = bt.badElecs;
        if isstruct(be)
            f2chk = {'noisyElecs','flatPSDElecs', ...
                     'badImpedanceElecs','declaredBadElectrodes'};
            for j = 1:numel(f2chk)
                if isfield(be, f2chk{j}) && ~isempty(be.(f2chk{j}))
                    badE = [badE; double(be.(f2chk{j})(:))]; %#ok<AGROW>
                end
            end
        end
    end
    S.badElectrodes = unique(badE(:))';

    % --------------------------------------------------------------- %
    %  3. impedance (optional)                                          %
    % --------------------------------------------------------------- %
    if exist(fullfile(folder,'impedanceData.mat'),'file')
        imp = load(fullfile(folder,'impedanceData.mat'));
        S.electrodeImpedances = double(imp.electrodeImpedances(:));
    else
        S.electrodeImpedances = [];
    end

    % --------------------------------------------------------------- %
    %  4. per-electrode data + topoplot positions                       %
    % --------------------------------------------------------------- %
    % Every electrode file's analogData is shaped [nTrials x nSamples].
    % We INSIST that the second dimension matches numel(S.timeVals);
    % otherwise the time-to-sample mapping would be silently wrong.
    nSamp = numel(S.timeVals);
    S.data          = NaN(Neeg, S.totalTrials, nSamp, 'single');
    S.elecAvailable = false(1, Neeg);
    S.elecXY        = NaN(Neeg, 2);

    for e = 1:Neeg
        f = fullfile(folder, sprintf('elec%d.mat', e));
        if exist(f,'file')
            d = load(f);
            % --- analogData ---
            if isfield(d,'analogData')
                ad = single(d.analogData);
                if size(ad,2) == nSamp
                    nTr = min(size(ad,1), S.totalTrials);
                    S.data(e, 1:nTr, :) = ad(1:nTr, :);
                    S.elecAvailable(e) = true;
                else
                    % Hard guarantee: timeVals MUST line up with samples.
                    warning('eegUtils:TimeMismatch', ...
                        ['elec%d.mat has %d samples but timeVals has %d - ' ...
                         'electrode skipped.'], e, size(ad,2), nSamp);
                end
            end
            % --- analogInfo (topoplot position) ---
            % Convention used here (matches EEGLAB):
            %     theta = azimuth in degrees (0 = front, +90 = right ear)
            %     radius = unit-disk radius (0 at vertex, ~0.5 at periphery)
            %     x = r * sin(theta),   y = r * cos(theta)
            if isfield(d,'analogInfo')
                ai = d.analogInfo;
                if isstruct(ai) && isfield(ai,'theta') && isfield(ai,'radius')
                    th = double(ai.theta);
                    r  = double(ai.radius);
                    if ~isempty(th) && ~isempty(r) && ~any(isnan([th(:); r(:)]))
                        S.elecXY(e,1) = r(1) * sind(th(1));
                        S.elecXY(e,2) = r(1) * cosd(th(1));
                    end
                end
            end
        end
    end

    % If any positions are still NaN (analogInfo missing or malformed) we
    % fall back to a hard-coded standard 10-10 map for the 64 actiCAP
    % labels.  This guarantees topoplots always render even when only a
    % subset of electrode files have analogInfo.
    if any(isnan(S.elecXY(:)))
        fb = standard1010Positions();
        for e = 1:Neeg
            if any(isnan(S.elecXY(e,:))) && isKey(fb, S.eegLabels{e})
                S.elecXY(e,:) = fb(S.eegLabels{e});
            end
        end
    end
end


%% ========================================================================
function [tIdx, tWin, tInfo] = getTimeWindowIdx(timeVals, window)
% GETTIMEWINDOWIDX  Map a time-window in seconds to sample indices.
% -------------------------------------------------------------------------
%   This is the single function every analysis uses to translate a request
%   like "the [-1 +1] s portion of each trial" into the actual sample
%   indices in the recorded data.  Doing it here once ensures that the
%   LFP, PAC and MI scripts ALL use the same timeVals-based selection.
%
%   INPUTS
%       timeVals : 1 x nSamp vector of times in seconds (e.g. -1.249..1.250)
%       window   : 1 x 2 vector [tStart tEnd] in seconds (inclusive)
%
%   OUTPUTS
%       tIdx     : 1 x nSamp logical mask, true for samples in [tStart,tEnd]
%       tWin     : the time vector restricted to the window (== timeVals(tIdx))
%       tInfo    : struct with diagnostics:
%                   .nSamples           number of samples kept
%                   .durationSec        actual duration kept (s)
%                   .actualStart/End    nearest sample times (s)
%                   .requestedStart/End the values you passed in
%                   .label              pretty-printed window string,
%                                        e.g. '[-1.000 to 1.000 s, 2001 samples]'
% -------------------------------------------------------------------------
    if ~isvector(timeVals)
        error('timeVals must be a 1-D vector.');
    end
    if numel(window) ~= 2 || window(2) < window(1)
        error('window must be [tStart tEnd] with tEnd >= tStart.');
    end

    timeVals = timeVals(:)';
    tIdx = (timeVals >= window(1)) & (timeVals <= window(2));
    tWin = timeVals(tIdx);

    if isempty(tWin)
        warning('eegUtils:EmptyWindow', ...
            'No samples in window [%.3f %.3f] s (timeVals span %.3f..%.3f s).', ...
            window(1), window(2), timeVals(1), timeVals(end));
    end

    tInfo.requestedStart = window(1);
    tInfo.requestedEnd   = window(2);
    tInfo.actualStart    = ifelse(isempty(tWin), NaN, tWin(1));
    tInfo.actualEnd      = ifelse(isempty(tWin), NaN, tWin(end));
    tInfo.nSamples       = numel(tWin);
    tInfo.durationSec    = ifelse(isempty(tWin), 0, tWin(end)-tWin(1));
    tInfo.label          = formatTimeWindow(tWin);
end


%% ========================================================================
function s = formatTimeWindow(tWin)
% FORMATTIMEWINDOW  Pretty-print an analysis window for plot titles.
%   Returns e.g. '[-1.000 to 1.000 s, 2001 samples]'
    if isempty(tWin)
        s = '[empty window]'; return;
    end
    s = sprintf('[%.3f to %.3f s, %d samples]', tWin(1), tWin(end), numel(tWin));
end


%% ========================================================================
function out = ifelse(cond, a, b)
% Tiny inline conditional (MATLAB has no built-in).
    if cond, out = a; else, out = b; end
end


%% ========================================================================
function fb = standard1010Positions()
% STANDARD1010POSITIONS  Hard-coded 10-10 (theta,radius) for 64 actiCAP labels.
% -------------------------------------------------------------------------
%   This is a fallback used by loadSubject when analogInfo.theta /
%   analogInfo.radius are missing.  Coordinates are converted to (x,y) on
%   a unit disk via x = r*sin(theta), y = r*cos(theta) - i.e. nose-up,
%   right-ear-positive-x.  The numbers below are the standard 10-10
%   positions used in EEGLAB's BESA template for these labels.
% -------------------------------------------------------------------------
labels = {'Fp1','Fz','F3','F7','FT9','FC5','FC1','C3','T7','TP9',...
          'CP5','CP1','Pz','P3','P7','O1','Oz','O2','P4','P8',...
          'TP10','CP6','CP2','Cz','C4','T8','FT10','FC6','FC2','F4',...
          'F8','Fp2','AF7','AF3','AFz','F1','F5','FT7','FC3','C1',...
          'C5','TP7','CP3','P1','P5','PO7','PO3','POz','PO4','PO8',...
          'P6','P2','CPz','CP4','TP8','C6','C2','FC4','FT8','F6',...
          'AF8','AF4','F2','Iz'};
tr = [
    -18 0.51; 0 0.25; -39 0.34; -54 0.51; -72 0.59; -69 0.51; -31 0.25; -90 0.33; -90 0.51; -108 0.59;
    -111 0.51; -149 0.25; 180 0.25; 141 0.34; 126 0.51; 162 0.51; 180 0.51; -162 0.51; -141 0.34; -126 0.51;
    108 0.59; 111 0.51; 149 0.25; 0 0; 90 0.33; 90 0.51; 72 0.59; 69 0.51; 31 0.25; 39 0.34;
    54 0.51; 18 0.51; -36 0.51; -21 0.34; 0 0.34; -23 0.20; -49 0.42; -65 0.51; -50 0.34; -90 0.17;
    -90 0.42; -90 0.59; -111 0.34; 159 0.20; 113 0.42; 144 0.51; 158 0.34; 180 0.34; -158 0.34; -144 0.51;
    -113 0.42; -159 0.20; 180 0.13; 111 0.34; 90 0.59; 90 0.42; 90 0.17; 50 0.34; 65 0.51; 49 0.42;
    36 0.51; 21 0.34; 23 0.20; 180 0.59];
fb = containers.Map();
for i = 1:numel(labels)
    th = tr(i,1); r = tr(i,2);
    fb(labels{i}) = [r*sind(th), r*cosd(th)];
end
end


%% ========================================================================
function [xy, ok] = getEEGPositions(S, elecIdx)
% GETEEGPOSITIONS  Look up topoplot (x,y) for a list of electrode indices.
xy = S.elecXY(elecIdx,:);
ok = ~any(isnan(xy),2);
end


%% ========================================================================
function y = bandpassFilter(x, Fs, lo, hi, order)
% BANDPASSFILTER  Zero-phase FIR bandpass.
% -------------------------------------------------------------------------
%   Uses fir1 (Hamming window) + filtfilt to avoid phase distortion.
%   Default order is Fs/2 taps (= 0.5 s for Fs=1000 Hz), which gives a
%   transition band of ~2 Hz at low frequencies.  filtfilt doubles the
%   effective order, so the realised attenuation is ~120 dB beyond the
%   transition band.
%
%   Edge handling: if hi >= Nyquist, it is clamped to Nyquist-1 Hz so the
%   filter design call doesn't error out.  Likewise lo <= 0 is clamped.
% -------------------------------------------------------------------------
if nargin < 5, order = round(Fs/2); end
if hi >= Fs/2, hi = Fs/2 - 1; end
if lo <= 0,    lo = 0.5;       end
b = fir1(order, [lo hi]/(Fs/2), 'bandpass');
isVec = isvector(x);
if isVec, x = x(:); end
y = zeros(size(x));
for k = 1:size(x,2)
    y(:,k) = filtfilt(b, 1, double(x(:,k)));
end
if isVec, y = reshape(y, size(x)); end
end


%% ========================================================================
function topoplotEEG(values, xy, varargin)
% TOPOPLOTEEG  Lightweight topographic plot WITHOUT the EEGLAB dependency.
% -------------------------------------------------------------------------
%   Renders a head outline, interpolates the per-electrode VALUES onto a
%   200x200 grid, masks to a circle of radius 0.5, and overlays the
%   electrode positions.  Designed to look the same across MATLAB
%   versions.
%
%   INPUTS
%       values : Neeg x 1 numeric (NaN-tolerant)
%       xy     : Neeg x 2 (x, y) on the unit disk (loadSubject's S.elecXY)
%
%   NAME-VALUE PAIRS
%       'labels'      cell array of channel labels (for showLabels=true)
%       'highlight'   electrode indices to draw with white-filled markers
%       'cmap'        colormap name (string) OR Nx3 matrix
%       'clim'        [cmin cmax] color limits
%       'title'       axis title
%       'showLabels'  print channel labels next to dots (default false)
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p,'labels',[]);
addParameter(p,'highlight',[]);
addParameter(p,'cmap','parula');
addParameter(p,'clim',[]);
addParameter(p,'title','');
addParameter(p,'showLabels',false);
parse(p,varargin{:});
opt = p.Results;

values = values(:);
keep = ~any(isnan(xy),2) & ~isnan(values);
x = xy(keep,1); y = xy(keep,2); v = values(keep);

if numel(v) < 4
    text(0,0,'Not enough valid electrodes for topoplot', ...
         'HorizontalAlignment','center');
    axis off; return;
end

% Interpolate on a regular grid clipped to the head circle (radius 0.5).
res = 200;
[xq,yq] = meshgrid(linspace(-0.55,0.55,res), linspace(-0.55,0.55,res));
F = scatteredInterpolant(x, y, double(v), 'natural', 'none');
Vq = F(xq, yq);
mask = (xq.^2 + yq.^2) <= 0.5^2;
Vq(~mask) = NaN;

imagesc(xq(1,:), yq(:,1), Vq); axis xy image off; hold on;
colormap(gca, opt.cmap);
if ~isempty(opt.clim), caxis(opt.clim); end

% Head outline + nose + ears.
th = linspace(0, 2*pi, 200);
plot(0.5*cos(th), 0.5*sin(th), 'k', 'LineWidth', 1.4);
plot([-0.05 0 0.05], [0.5 0.55 0.5], 'k', 'LineWidth', 1.4);
plot([-0.5 -0.55 -0.55 -0.5], [0.1 0.05 -0.05 -0.1], 'k', 'LineWidth', 1.2);
plot([ 0.5  0.55  0.55  0.5], [0.1 0.05 -0.05 -0.1], 'k', 'LineWidth', 1.2);

plot(x, y, 'k.', 'MarkerSize', 6);
if ~isempty(opt.highlight)
    keepMap = find(keep);
    [~, idxInKept] = ismember(opt.highlight, keepMap);
    idxInKept = idxInKept(idxInKept>0);
    if ~isempty(idxInKept)
        plot(x(idxInKept), y(idxInKept), 'ko', ...
             'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.2);
    end
end
if opt.showLabels && ~isempty(opt.labels)
    lbl = opt.labels(keep);
    for i = 1:numel(lbl)
        text(x(i), y(i)+0.03, lbl{i}, ...
             'HorizontalAlignment','center','FontSize',7);
    end
end
if ~isempty(opt.title), title(opt.title, 'FontWeight','bold'); end
end


%% ========================================================================
function saveFigDual(fig, baseName)
% SAVEFIGDUAL  Save a figure as both PNG (publication) and FIG (editable).
%   This pair lets you embed the PNG in LaTeX immediately AND still
%   re-tweak any plot detail later by opening the .fig in MATLAB.
set(fig,'PaperPositionMode','auto');
print(fig, [baseName '.png'], '-dpng', '-r300');
savefig(fig, [baseName '.fig']);
end


%% ========================================================================
function pAdj = fdrBH(p, ~)
% FDRBH  Benjamini-Hochberg FDR-corrected p-values.
% -------------------------------------------------------------------------
%   Standard step-up procedure that controls the false discovery rate at
%   level alpha (compare pAdj < alpha downstream).  Equivalent
%   to mafdr(p,'BHFDR',true) but works without the Statistics Toolbox.
% -------------------------------------------------------------------------
p = p(:);
[ps, ix] = sort(p);
m = numel(p);
adj = ps .* m ./ (1:m)';
adj = min(1, cummin(adj, 'reverse'));   % monotone non-decreasing in rank
pAdj = nan(size(p));
pAdj(ix) = adj;
end


%% ========================================================================
function elecs = defaultFrontoParietal()
% DEFAULTFRONTOPARIETAL  10 fronto-parietal electrodes (FPN-style).
% -------------------------------------------------------------------------
%   Returned in the order [frontal | parietal] so plots that show them
%   left-to-right give an anatomically intuitive layout.
%
%   Selection rationale
%   -------------------
%   The fronto-parietal network (FPN) underlies attentional control and
%   working memory.  In EEG it is conventionally probed at the dorsal
%   midline-to-lateral frontal sites and the parietal sites overlying
%   posterior parietal cortex.  We pick five of each:
%
%       Frontal :  F3 (3),  Fz (2),  F2 (63),  AFz (35),  AF4 (62)
%       Parietal:  P3 (14), Pz (13), P4 (19),  P1  (44),  P2  (52)
%
%   IMPORTANT - bad-electrode awareness
%   -----------------------------------
%   In the sample subject's badTrials_wo_v8.mat, the union of
%   noisyElecs ([5 10 29 30]) and flatPSDElecs ([34]) is
%       [5 (FT9), 10 (TP9), 29 (FC2), 30 (F4), 34 (AF3)]
%   We therefore deliberately AVOID indices 5, 10, 29, 30 and 34 in this
%   default set.  The substitutions are:
%       F4 (30)  ->  F2 (63)    [right-frontal, slightly more medial]
%       AF3 (34) ->  AFz (35)   [midline anterior-frontal]
%   so the canonical set still spans bilateral dorsolateral frontal
%   cortex and posterior parietal cortex without forcing the loader to
%   discard any focus electrode in this subject.
%
%   For OTHER subjects, different electrodes may be flagged bad.  The
%   downstream scripts handle that automatically by intersecting this
%   canonical set with each subject's good-electrode list:
%       focus_for_subject = intersect(defaultFrontoParietal, goodElecs)
%   so any subject-specific bad electrode is dropped at run-time.
% -------------------------------------------------------------------------
elecs = [3 2 63 35 62  14 13 19 44 52];   %#ok<NBRAK> kept on one line
end


%% ========================================================================
function printDivider(msg)
% PRINTDIVIDER  Cosmetic CLI banner, used by every analysis on entry/exit.
fprintf('\n========================================================\n');
fprintf('  %s\n', msg);
fprintf('========================================================\n');
end
