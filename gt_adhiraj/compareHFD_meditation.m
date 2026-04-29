% =========================================================================
% compareHFD_meditation.m
% Compare meditation-condition ΔHFD between meditator and control.
% =========================================================================

% ================== CONFIG ==================

baseDir = 'C:\Users\lenovo\OneDrive - Indian Institute of Science\Documents\qtbio_MATLAB\courses\NSP\code\provided\Project_Master';

subjects = {'013AR','064PK'};
subjectLabels = {'013AR meditator','064PK control'};

resultDirs = {
    fullfile(baseDir, 'results_013AR')
    fullfile(baseDir, 'results_064PK')
};

medConds = {'M1','M2'};

HP_order = {'O1','Oz','O2','PO3','POz','PO4','P1','P3','P2'};

% subject × meditation condition
deltaMed = zeros(numel(subjects), numel(medConds));
pMed     = zeros(numel(subjects), numel(medConds));

% ================== LOAD + COMPUTE ==================

for s = 1:numel(subjects)

    subjectID = subjects{s};

    for m = 1:numel(medConds)

        cond = medConds{m};

        fileName = ['HFD_' subjectID '_' cond '_ep_v8.mat'];
        fullPath = fullfile(resultDirs{s}, fileName);

        if ~exist(fullPath, 'file')
            error('Missing file: %s', fullPath);
        end

        load(fullPath, 'HFD_BL', 'HFD_ST', 'all_labels');

        hp_idx = find(ismember(all_labels, HP_order));

        % Pool across HP electrodes and trials
        hp_BL = reshape(HFD_BL(hp_idx,:), 1, []);
        hp_ST = reshape(HFD_ST(hp_idx,:), 1, []);

        deltaVals = hp_ST - hp_BL;

        deltaMed(s,m) = median(deltaVals);
        pMed(s,m) = signrank(hp_BL', hp_ST');

    end
end

% ================== PRINT TABLES ==================

fprintf('\nMedian HP pooled ΔHFD during meditation conditions:\n');
T_delta = array2table(deltaMed, ...
    'VariableNames', medConds, ...
    'RowNames', subjects);
disp(T_delta);

fprintf('p-values:\n');
T_p = array2table(pMed, ...
    'VariableNames', medConds, ...
    'RowNames', subjects);
disp(T_p);

% ================== PLOT ==================

fig = figure('Name','Meditation HFD comparison', ...
             'Position',[100 100 800 500]);
hold on;

x = 1:numel(medConds);

plot(x, deltaMed(1,:), 'o-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'DisplayName', subjectLabels{1});

plot(x, deltaMed(2,:), 's-', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 8, ...
    'DisplayName', subjectLabels{2});

yline(0, 'k--', 'LineWidth', 1);

set(gca, ...
    'XTick', x, ...
    'XTickLabel', medConds, ...
    'FontSize', 11);

xlabel('Meditation condition');
ylabel('Median HP pooled ΔHFD (ST − BL)');
title('Meditation-related EEG complexity modulation');

legend('Location','best');
box off;

% ================== SAVE ==================

compareDir = fullfile(baseDir, 'results_comparisons');

if ~exist(compareDir, 'dir')
    mkdir(compareDir);
end

% Save figures
outPNG = fullfile(compareDir, 'Compare_HFD_Meditation_013AR_vs_064PK.png');
outFIG = fullfile(compareDir, 'Compare_HFD_Meditation_013AR_vs_064PK.fig');

saveas(fig, outPNG);
saveas(fig, outFIG);

fprintf('\nSaved meditation comparison figure:\n%s\n%s\n', outPNG, outFIG);

% Save readable CSV tables
csvDelta = fullfile(compareDir, 'Compare_HFD_Meditation_delta.csv');
csvPvals = fullfile(compareDir, 'Compare_HFD_Meditation_pvalues.csv');

writetable(T_delta, csvDelta, 'WriteRowNames', true);
writetable(T_p, csvPvals, 'WriteRowNames', true);

fprintf('Saved CSV tables:\n%s\n%s\n', csvDelta, csvPvals);