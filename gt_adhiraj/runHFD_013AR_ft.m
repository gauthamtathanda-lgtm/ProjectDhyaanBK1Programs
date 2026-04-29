% =========================================================================
% runHFD_013AR.m  —  UPDATED for FieldTrip-preprocessed data (013AR_ft)
% G1 Project: Fractal Analysis of EEG in Meditators
%
% DATA FORMAT (013AR_ft, *_ep_v8.mat):
%   Each file contains a FieldTrip epoched struct with fields:
%     data.trial  — {1 × nTrials} cell, each cell = (64 × 2500) double
%     data.time   — {1 × nTrials} cell, each cell = time vector (s)
%     data.label  — {64 × 1} cell of channel name strings
%     data.fsample — scalar (1000 Hz)
%     data.timeVals — (1 × 2500) shared time vector
%     data.badElecs — electrode indices flagged as bad (1-indexed)
%     data.numTrials / numGoodTrials — scalars
%
% KEY DIFFERENCES FROM OLD PIPELINE:
%   - Old: individual elecN.mat files at 2000 Hz, 120 trials, (trials×2500)
%   - New: single FieldTrip struct per condition at 1000 Hz, variable trials
%   - New: all 64 channels in one file; bad electrodes listed in data.badElecs
%   - New: no separate filtering needed — data already preprocessed (v8)
%   - New: BL/ST windows use 500 samples each (same 500ms, adjusted for 1kHz)
%   - New: kmax=10 broadband (unchanged), kmax=5 for short windows (unchanged)
%   - New: 87 G1 trials → 3 spf × 29 trials per spf (not 40/40/40)
%   - New: Oz cross-check done on first good trial of G1 baseline
%
% BEFORE RUNNING:
%   1. Place all *_ep_v8.mat files in dataDir (or update dataDir below)
%   2. Ensure computeHFD.m and computeDFA.m are on the MATLAB path
%   3. Chronux is only needed for the PSD section (Step 9)
% =========================================================================


% ================== SUBJECT CONFIG ==================

subjectID = '064PK';   % use '013AR' for meditator, '064PK' for control

baseDir = 'C:\Users\lenovo\OneDrive - Indian Institute of Science\Documents\qtbio_MATLAB\courses\NSP\code\provided';

dataDir   = fullfile(baseDir, [subjectID '_ft']);
outputDir = fullfile(baseDir, 'Project_Master', ['results_' subjectID]);
codeDir   = fileparts(mfilename('fullpath'));

addpath(genpath('/path/to/chronux_2_12'));
addpath(codeDir);

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

conditionFiles = {
    'EC1_ep_v8.mat'
    'EC2_ep_v8.mat'
    'EO1_ep_v8.mat'
    'EO2_ep_v8.mat'
    'G1_ep_v8.mat'
    'G2_ep_v8.mat'
    'M1_ep_v8.mat'
    'M2_ep_v8.mat'
};

fprintf('Subject: %s\n', subjectID);
fprintf('Data folder: %s\n', dataDir);
fprintf('Output folder: %s\n', outputDir);
fprintf('Data folder found. Ready.\n');

% ================== LOOP START ==================

for c = 1:numel(conditionFiles)

    gFile = fullfile(dataDir, conditionFiles{c});
    [~,name,~] = fileparts(gFile);

    fprintf('\n\n====================================================\n');
    fprintf('Running condition %d/%d: %s\n', c, numel(conditionFiles), conditionFiles{c});
    fprintf('====================================================\n');

    if ~exist(gFile,'file')
        error('%s not found in dataDir:\n  %s\nCheck your path.', conditionFiles{c}, dataDir);
    end

% ── 1. PARAMETERS ────────────────────────────────────────────────────────

Fs      = 1000;          % sampling rate (Hz) — new data is 1 kHz (was 2 kHz)
nSamp   = 2500;          % samples per trial (epoch: -1.249 to +1.250 s)
kmax_broad = 10;         % kmax for broadband HFD (1–90 Hz)
kmax_short = 5;          % kmax for 100ms windows

% Epoch windows — defined after data load using ftData.timeVals
% (avoids floating-point edge issue from linspace that caused 499 not 500 samples)

fprintf('Fs = %d Hz\n', Fs);

% ── 2. ELECTRODE MAPPING ─────────────────────────────────────────────────
% HP (high-priority) occipito-parietal group — Murty et al. 2020
% Channel names identical across all *_ep_v8.mat files (64-ch BrainAmp layout)
HP_labels   = {'O1','Oz','O2','PO3','POz','PO4','P1','P3','P2'};
NONHP_labels= {'Cz','C4','Fz','F3'};   % frontal/central comparators
% Note: old NONHP were C2, FC4, FT8, Fp1 — those are now marked as badElecs
% in 013AR_ft (FC4→ch57, C2→ch58 flagged). Using clean alternatives instead.

all_labels = [HP_labels, NONHP_labels];
nAll       = numel(all_labels);
nHP        = numel(HP_labels);

% ── 3. LOAD G1 DATA ──────────────────────────────────────────────────────
fprintf('\nLoading %s ...\n', [name '.mat']);
raw     = load(gFile, 'data', 'numGoodTrials');
ftData  = raw.data;           % FieldTrip struct
nTrials = double(raw.numGoodTrials);

% Build timeVec from the actual stored timeVals — exact, no rounding error
timeVec = double(ftData.timeVals(:)');
bl_idx  = timeVec >= -0.5  & timeVec <= -0.001;  % baseline: -500 to 0 ms  (500 samp)
st_idx  = timeVec >=  0.25 & timeVec <  0.750;  % stimulus: 250 to 750 ms (500 samp)
fprintf('BL window: %d samples (%.0f ms)\n', sum(bl_idx), sum(bl_idx)/Fs*1000);
fprintf('ST window: %d samples (%.0f ms)\n', sum(st_idx), sum(st_idx)/Fs*1000);

% Extract channel label list from the FieldTrip struct
ft_labels = cellfun(@(x) char(x), ftData.label, 'UniformOutput', false);
ft_labels = ft_labels(:);     % ensure column cell array

% Bad electrode names (1-indexed in data.badElecs)
badIdx      = ftData.badElecs(:);
badNames    = ft_labels(badIdx);
fprintf('numGoodTrials = %d\n', nTrials);
fprintf('Bad electrodes (%d): %s\n', numel(badNames), strjoin(badNames, ', '));

% Warn if any HP or NONHP electrodes are bad
for i = 1:nAll
    if any(strcmp(badNames, all_labels{i}))
        warning('Electrode %s is flagged as bad — results may be unreliable.', all_labels{i});
    end
end

% Map our electrode labels to FieldTrip channel indices
elec_idx = zeros(1, nAll);
for i = 1:nAll
    idx = find(strcmp(ft_labels, all_labels{i}));
    if isempty(idx)
        error('Electrode %s not found in data.label.', all_labels{i});
    end
    elec_idx(i) = idx;
end

% ── 4. BUILD DATA MATRIX (channels × trials × samples) ──────────────────
% FieldTrip trial cell: each ftData.trial{tr} = (64 × 2500)
% We extract only the electrodes we need → (nAll × nTrials × nSamp)
fprintf('\nExtracting electrode data ...\n');
eegData = zeros(nAll, nTrials, nSamp);
for tr = 1:nTrials
    mat_tr = ftData.trial{tr};          % (64 × 2500) double
    for ei = 1:nAll
        eegData(ei, tr, :) = mat_tr(elec_idx(ei), :);
    end
end

% ── 5. FILTERS ───────────────────────────────────────────────────────────
% Data is already preprocessed, but we re-apply the same bandpass + notch
% for consistency with the analysis framework (Aggarwal & Ray 2025).
% Filters designed at NEW Fs = 1000 Hz (was 2000 Hz — cutoff ratios identical)
fprintf('Designing filters at Fs=%d Hz ...\n', Fs);

[b_bp, a_bp] = cheby1(4, 0.2, [1 90]  / (Fs/2), 'bandpass');  % 1–90 Hz
[b_bs, a_bs] = cheby1(4, 0.2, [48 52] / (Fs/2), 'stop');       % 50 Hz notch

% Apply filters
eegFilt = zeros(size(eegData));
for ei = 1:nAll
    for tr = 1:nTrials
        s = filtfilt(b_bs, a_bs, squeeze(eegData(ei,tr,:)));
        eegFilt(ei, tr, :) = filtfilt(b_bp, a_bp, s);
    end
end
fprintf('Filtering done.\n');

% ── 6. COMPUTE HFD ───────────────────────────────────────────────────────
fprintf('\nComputing HFD (kmax=%d broadband)...\n', kmax_broad);

HFD_BL = zeros(nAll, nTrials);
HFD_ST = zeros(nAll, nTrials);

for ei = 1:nAll
    for tr = 1:nTrials
        sig = squeeze(eegFilt(ei, tr, :));
        HFD_BL(ei, tr) = computeHFD(sig(bl_idx), kmax_broad);
        HFD_ST(ei, tr) = computeHFD(sig(st_idx), kmax_broad);
    end
    dHFD_med = median(HFD_ST(ei,:) - HFD_BL(ei,:));
    p = signrank(HFD_BL(ei,:)', HFD_ST(ei,:)');
    fprintf('  %-5s: BL=%.4f  ST=%.4f  Δ=%+.4f  p=%.2e\n', ...
            all_labels{ei}, median(HFD_BL(ei,:)), median(HFD_ST(ei,:)), dHFD_med, p);
end

dHFD = HFD_ST - HFD_BL;

% ── 7. HFD CROSS-CHECK: HFD + H ≈ 1.98 ─────────────────────────────────
% Validates kmax selection — Aggarwal & Ray (2025) Section 3.1.2.1
fprintf('\nCross-check: HFD + H (expect ≈ 1.98)...\n');
oz_pos = find(strcmp(all_labels,'Oz'));
H_vals = zeros(1, nTrials);
for tr = 1:nTrials
    % Use full filtered trial (2500 samples) for DFA — Aggarwal & Ray (2025)
    % use the full epoch signal for H estimation; 500-sample BL window is
    % too short for reliable DFA (nMax = 125 gives only ~9 scale points).
    H_vals(tr) = computeDFA(squeeze(eegFilt(oz_pos, tr, :)));
end
crosscheck = mean(HFD_BL(oz_pos,:)) + mean(H_vals);
fprintf('  HFD + H = %.4f  (expect ~1.98)\n', crosscheck);

% ── 8. HP GROUP STATISTICS ───────────────────────────────────────────────
hp_idx     = 1:nHP;    % HP electrodes are first in all_labels
hp_BL_flat = reshape(HFD_BL(hp_idx,:), 1, []);
hp_ST_flat = reshape(HFD_ST(hp_idx,:), 1, []);
p_hp       = signrank(hp_BL_flat', hp_ST_flat');

fprintf('\n=== HP GROUP POOLED ===\n');
fprintf('  BL   = %.4f ± %.4f\n', median(hp_BL_flat), std(hp_BL_flat));
fprintf('  ST   = %.4f ± %.4f\n', median(hp_ST_flat), std(hp_ST_flat));
fprintf('  ΔHFD = %+.4f\n', median(hp_ST_flat - hp_BL_flat));
fprintf('  p    = %.2e\n', p_hp);

% ── 9. ΔHFD BY SPATIAL FREQUENCY ────────────────────────────────────────
% G1 has 87 trials = 3 spf × 29 trials per spf
% (NEW: 29 per spf, not 40/40/40 as in old 120-trial dataset)
% Assumed ordering: trials 1–29 = 1cpd, 30–58 = 2cpd, 59–87 = 4cpd
%   *** Verify this ordering with your paradigm log if available ***
nPerSpf  = floor(nTrials / 3);          % = 29
spf_vals = [ones(1,nPerSpf), 2*ones(1,nPerSpf), 4*ones(1,nPerSpf)];
spf_list = [1 2 4];
nUsed    = 3 * nPerSpf;                 % 87 (all trials used here)

fprintf('\n=== ΔHFD BY SPATIAL FREQUENCY (HP pooled) ===\n');
fprintf('    (nPerSpf = %d, total used = %d / %d)\n', nPerSpf, nUsed, nTrials);
spf_meds = zeros(1,3);
for si = 1:3
    idx  = spf_vals == spf_list(si);
    vals = reshape(dHFD(hp_idx, idx), 1, []);
    bl_v = reshape(HFD_BL(hp_idx, idx), 1, []);
    st_v = reshape(HFD_ST(hp_idx, idx), 1, []);
    p_spf = signrank(bl_v', st_v');
    spf_meds(si) = median(vals);
    fprintf('  %d cpd: ΔHFD=%+.5f  p=%.2e\n', spf_list(si), spf_meds(si), p_spf);
end

% ── 10. TIME-RESOLVED HFD (100ms windows) ────────────────────────────────
fprintf('\nComputing time-resolved HFD (100ms windows, kmax=%d)...\n', kmax_short);
win_len = round(0.1 * Fs);              % 100 samples at 1 kHz = 100ms
nWins   = floor(nSamp / win_len);       % 25 windows across 2500ms epoch
tc_time = linspace(timeVec(1) + 0.05, timeVec(end) - 0.05, nWins);

tc_all = zeros(nHP, nTrials, nWins);
for i = 1:nHP
    for tr = 1:nTrials
        sig = squeeze(eegFilt(i, tr, :));
        for w = 1:nWins
            seg = sig((w-1)*win_len+1 : w*win_len);
            tc_all(i, tr, w) = computeHFD(seg, kmax_short);
        end
    end
    fprintf('  %s done\n', HP_labels{i});
end

tc_hp  = squeeze(mean(tc_all, 1));   % (nTrials × nWins) — avg over HP elecs
tc_med = median(tc_hp, 1);
tc_sem = std(tc_hp, [], 1) / sqrt(nTrials);

% ── 11. POWER SPECTRUM (Oz electrode) ───────────────────────────────────
% Requires Chronux (mtspectrumc). Skip this section if not available.
if exist('mtspectrumc','file')
    oz_data = squeeze(eegFilt(oz_pos, :, :));    % (nTrials × nSamp)
    params.tapers = [3 5]; params.Fs = Fs;
    params.fpass  = [1 90]; params.pad = 1;

    [S_bl, f] = mtspectrumc(oz_data(:, bl_idx)', params);
    [S_st, ~] = mtspectrumc(oz_data(:, st_idx)', params);
    psd_bl = mean(S_bl, 2);
    psd_st = mean(S_st, 2);
    fprintf('\nPSD computed (Chronux).\n');
else
    warning('Chronux not found — using pwelch for PSD instead.');
    oz_data = squeeze(eegFilt(oz_pos, :, :));
    f    = 1:90;
    psd_bl = zeros(numel(f),1);
    psd_st = zeros(numel(f),1);
    for tr = 1:nTrials
        [p1,fp] = pwelch(oz_data(tr, bl_idx), 256, 128, f, Fs);
        [p2,~]  = pwelch(oz_data(tr, st_idx), 256, 128, f, Fs);
        psd_bl  = psd_bl + p1(:);
        psd_st  = psd_st + p2(:);
    end
    psd_bl = psd_bl / nTrials;
    psd_st = psd_st / nTrials;
    f = f(:);
    fprintf('\nPSD computed (pwelch fallback).\n');
end

% ── 12. SAVE RESULTS ─────────────────────────────────────────────────────
% Save one result file per condition.

outFile = fullfile(outputDir, ['HFD_' subjectID '_' name '.mat']);

save(outFile, ...
     'HFD_BL','HFD_ST','dHFD','all_labels','HP_labels', ...
     'hp_BL_flat','hp_ST_flat','p_hp', ...
     'tc_med','tc_sem','tc_time', ...
     'psd_bl','psd_st','f', ...
     'spf_meds','timeVec','Fs', ...
     'nTrials','nPerSpf','spf_vals', ...
     'eegFilt', '-v7.3');

fprintf('\nResults saved to:\n  %s\n', outFile);
fprintf('Finished condition: %s\n', conditionFiles{c});

end