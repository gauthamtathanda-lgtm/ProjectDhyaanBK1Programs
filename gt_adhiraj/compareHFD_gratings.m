% =========================================================================
% compareHFD_gratings.m
% Compare spatial-frequency-dependent ΔHFD for G1/G2 in meditator vs control.
% =========================================================================

% ================== CONFIG ==================

baseDir = 'C:\Users\lenovo\OneDrive - Indian Institute of Science\Documents\qtbio_MATLAB\courses\NSP\code\provided\Project_Master';

subjects = {'013AR','064PK'};
subjectLabels = {'013AR meditator','064PK control'};

resultDirs = {
    fullfile(baseDir, 'results_013AR')
    fullfile(baseDir, 'results_064PK')
};

gratingConds = {'G1','G2'};
spfLabels = {'1 cpd','2 cpd','4 cpd'};
spfVals = [1 2 4];

HP_order = {'O1','Oz','O2','PO3','POz','PO4','P1','P3','P2'};

% dimensions: subject × grating condition × spatial frequency
deltaSPF = zeros(numel(subjects), numel(gratingConds), numel(spfVals));
pSPF     = zeros(numel(subjects), numel(gratingConds), numel(spfVals));

% ================== LOAD + COMPUTE ==================

for s = 1:numel(subjects)

    subjectID = subjects{s};

    for g = 1:numel(gratingConds)

        cond = gratingConds{g};
        fileName = ['HFD_' subjectID '_' cond '_ep_v8.mat'];
        fullPath = fullfile(resultDirs{s}, fileName);

        if ~exist(fullPath, 'file')
            error('Missing file: %s', fullPath);
        end

        load(fullPath, 'HFD_BL', 'HFD_ST', 'dHFD', ...
             'all_labels', 'nPerSpf');

        hp_idx = find(ismember(all_labels, HP_order));

        spfTrialVals = [ones(1,nPerSpf), ...
                        2*ones(1,nPerSpf), ...
                        4*ones(1,nPerSpf)];

        nUsed = 3 * nPerSpf;

        for si = 1:numel(spfVals)

            idx = spfTrialVals == spfVals(si);

            vals = reshape(dHFD(hp_idx, idx), 1, []);
            bl_v = reshape(HFD_BL(hp_idx, idx), 1, []);
            st_v = reshape(HFD_ST(hp_idx, idx), 1, []);

            deltaSPF(s,g,si) = median(vals);
            pSPF(s,g,si) = signrank(bl_v', st_v');

        end

        fprintf('%s %s: used %d/%d trials for SPF analysis\n', ...
                subjectID, cond, nUsed, size(HFD_BL,2));

    end
end

% ================== PRINT TABLES ==================

for g = 1:numel(gratingConds)

    fprintf('\nMedian HP pooled ΔHFD by SPF — %s\n', gratingConds{g});
    T_delta = array2table(squeeze(deltaSPF(:,g,:)), ...
        'VariableNames', spfLabels, ...
        'RowNames', subjects);
    disp(T_delta);

    fprintf('p-values by SPF — %s\n', gratingConds{g});
    T_p = array2table(squeeze(pSPF(:,g,:)), ...
        'VariableNames', spfLabels, ...
        'RowNames', subjects);
    disp(T_p);

end

% ================== PLOT ==================

fig = figure('Name','Grating SPF comparison', ...
             'Position',[100 100 1100 450]);

for g = 1:numel(gratingConds)

    subplot(1,2,g); hold on;

    x = 1:numel(spfVals);

    plot(x, squeeze(deltaSPF(1,g,:)), 'o-', ...
        'LineWidth', 2.5, ...
        'MarkerSize', 8, ...
        'DisplayName', subjectLabels{1});

    plot(x, squeeze(deltaSPF(2,g,:)), 's-', ...
        'LineWidth', 2.5, ...
        'MarkerSize', 8, ...
        'DisplayName', subjectLabels{2});

    yline(0, 'k--', 'LineWidth', 1);

    set(gca, ...
        'XTick', x, ...
        'XTickLabel', spfLabels, ...
        'FontSize', 11);

    xlabel('Spatial frequency');
    ylabel('Median HP pooled ΔHFD');
    title(sprintf('%s: spatial-frequency modulation', gratingConds{g}));

    legend('Location','best');
    box off;

end

sgtitle('Grating-driven HFD modulation: meditator vs control', ...
        'FontSize', 13, ...
        'FontWeight', 'bold');

% ================== SAVE ==================

compareDir = fullfile(baseDir, 'results_comparisons');

if ~exist(compareDir, 'dir')
    mkdir(compareDir);
end

outPNG = fullfile(compareDir, 'Compare_HFD_Gratings_SPF_013AR_vs_064PK.png');
outFIG = fullfile(compareDir, 'Compare_HFD_Gratings_SPF_013AR_vs_064PK.fig');

saveas(fig, outPNG);
saveas(fig, outFIG);

fprintf('\nSaved grating comparison figure:\n%s\n%s\n', outPNG, outFIG);

% ---- save readable CSV tables ----
for g = 1:numel(gratingConds)

    T_delta = array2table(squeeze(deltaSPF(:,g,:)), ...
        'VariableNames', spfLabels, ...
        'RowNames', subjects);

    T_p = array2table(squeeze(pSPF(:,g,:)), ...
        'VariableNames', spfLabels, ...
        'RowNames', subjects);

    csvDelta = fullfile(compareDir, ...
        ['Compare_HFD_Gratings_' gratingConds{g} '_delta.csv']);

    csvPvals = fullfile(compareDir, ...
        ['Compare_HFD_Gratings_' gratingConds{g} '_pvalues.csv']);

    writetable(T_delta, csvDelta, 'WriteRowNames', true);
    writetable(T_p, csvPvals, 'WriteRowNames', true);

    fprintf('Saved CSV tables:\n%s\n%s\n', csvDelta, csvPvals);

end