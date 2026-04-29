% =========================================================================
% compareHFD_subjects.m
% Compare HP pooled ΔHFD across conditions for meditator vs control.
% =========================================================================

% ================== CONFIG ==================

baseDir = 'C:\Users\lenovo\OneDrive - Indian Institute of Science\Documents\qtbio_MATLAB\courses\NSP\code\provided\Project_Master';

subjects = {'013AR','064PK'};

resultDirs = {
    fullfile(baseDir, 'results_013AR')
    fullfile(baseDir, 'results_064PK')
};

conditions = {'EC1','EC2','EO1','EO2','G1','G2','M1','M2'};

% ================== STORAGE ==================

deltaMat = zeros(numel(subjects), numel(conditions));
pMat     = zeros(numel(subjects), numel(conditions));

% ================== LOAD + COMPUTE ==================

for s = 1:numel(subjects)

    subjectID = subjects{s};

    for c = 1:numel(conditions)

        cond = conditions{c};

        fileName = ['HFD_' subjectID '_' cond '_ep_v8.mat'];
        fullPath = fullfile(resultDirs{s}, fileName);

        if ~exist(fullPath, 'file')
            error('Missing file: %s', fullPath);
        end

        load(fullPath, 'HFD_BL', 'HFD_ST', 'all_labels');

        % HP electrodes: occipito-parietal group
        HP_order = {'O1','Oz','O2','PO3','POz','PO4','P1','P3','P2'};
        hp_idx = find(ismember(all_labels, HP_order));

        % Pool across HP electrodes and trials
        hp_BL = reshape(HFD_BL(hp_idx,:), 1, []);
        hp_ST = reshape(HFD_ST(hp_idx,:), 1, []);

        deltaVals = hp_ST - hp_BL;

        % Store median ΔHFD and paired sign-rank p-value
        deltaMat(s,c) = median(deltaVals);
        pMat(s,c) = signrank(hp_BL', hp_ST');

    end
end

% ================== PRINT TABLES ==================

disp('Median HP pooled ΔHFD:');
T_delta = array2table(deltaMat, ...
    'VariableNames', conditions, ...
    'RowNames', subjects);
disp(T_delta);

disp('p-values:');
T_p = array2table(pMat, ...
    'VariableNames', conditions, ...
    'RowNames', subjects);
disp(T_p);

% ================== PLOT ==================

fig = figure('Name','HFD comparison: 013AR vs 064PK', ...
             'Position',[100 100 1000 500]);
hold on;

x = 1:numel(conditions);

plot(x, deltaMat(1,:), 'o-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'DisplayName', '013AR meditator');

plot(x, deltaMat(2,:), 's-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'DisplayName', '064PK control');

yline(0, 'k--', 'LineWidth', 1);

set(gca, ...
    'XTick', x, ...
    'XTickLabel', conditions, ...
    'FontSize', 11);

xlabel('Condition');
ylabel('Median HP pooled ΔHFD (Stimulus - Baseline)');
title('EEG complexity modulation: meditator vs control');

legend('Location','best');
box off;

% ================== SAVE ==================

compareDir = fullfile(baseDir, 'results_comparisons');

if ~exist(compareDir, 'dir')
    mkdir(compareDir);
end

% ---- save figures ----
outPNG = fullfile(compareDir, 'Compare_HFD_013AR_vs_064PK.png');
outFIG = fullfile(compareDir, 'Compare_HFD_013AR_vs_064PK.fig');

saveas(fig, outPNG);
saveas(fig, outFIG);

fprintf('\nSaved comparison figure:\n%s\n%s\n', outPNG, outFIG);

% ---- save readable CSV tables ----
T_delta = array2table(deltaMat, ...
    'VariableNames', conditions, ...
    'RowNames', subjects);

T_p = array2table(pMat, ...
    'VariableNames', conditions, ...
    'RowNames', subjects);

csvDelta = fullfile(compareDir, 'Compare_HFD_delta.csv');
csvPvals = fullfile(compareDir, 'Compare_HFD_pvalues.csv');

writetable(T_delta, csvDelta, 'WriteRowNames', true);
writetable(T_p, csvPvals, 'WriteRowNames', true);

fprintf('Saved CSV tables:\n%s\n%s\n', csvDelta, csvPvals);