% =========================================================================
% plotHFD_013AR_ft.m
% Batch plot generator for FieldTrip-preprocessed 013AR_ft HFD results.
% Run AFTER runHFD_013AR_ft.m has generated all HFD_<subjectID>_*.mat files.
% =========================================================================

% ================== CONFIG ==================

subjectID = '064PK';

baseDir = 'C:\Users\lenovo\OneDrive - Indian Institute of Science\Documents\qtbio_MATLAB\courses\NSP\code\provided';

outputDir = fullfile(baseDir, 'Project_Master', ['results_' subjectID]);

resultFiles = {
    ['HFD_' subjectID '_EC1_ep_v8.mat']
    ['HFD_' subjectID '_EC2_ep_v8.mat']
    ['HFD_' subjectID '_EO1_ep_v8.mat']
    ['HFD_' subjectID '_EO2_ep_v8.mat']
    ['HFD_' subjectID '_G1_ep_v8.mat']
    ['HFD_' subjectID '_G2_ep_v8.mat']
    ['HFD_' subjectID '_M1_ep_v8.mat']
    ['HFD_' subjectID '_M2_ep_v8.mat']
};

fprintf('Starting batch plotting...\n');

for r = 1:numel(resultFiles)

    resultFile = resultFiles{r};
    fullPath = fullfile(outputDir, resultFile);

    if ~exist(fullPath, 'file')
        warning('%s not found. Skipping.', resultFile);
        continue;
    end

    tag = erase(resultFile, ['HFD_' subjectID '_']);
    tag = erase(tag, '.mat');

    fprintf('\n========================================\n');
    fprintf('Plotting %d/%d: %s\n', r, numel(resultFiles), resultFile);
    fprintf('========================================\n');

    load(fullPath);

    % ── Recompute per-electrode p-values ─────────────────────────────────
    nElec = size(HFD_BL, 1);
    pvals = zeros(nElec, 1);

    for i = 1:nElec
        pvals(i) = signrank(HFD_BL(i,:)', HFD_ST(i,:)');
    end

    % HP group pooled
    HP_order   = {'O1','Oz','O2','PO3','POz','PO4','P1','P3','P2'};
    hp_idx_arr = find(ismember(all_labels, HP_order));

    hp_BL_flat = reshape(HFD_BL(hp_idx_arr,:), 1, []);
    hp_ST_flat = reshape(HFD_ST(hp_idx_arr,:), 1, []);
    p_hp = signrank(hp_BL_flat', hp_ST_flat');

    % Rebuild spatial-frequency vector safely
    spf_list = [1 2 4];
    nUsed = 3 * nPerSpf;
    spf_vals = [ones(1,nPerSpf), 2*ones(1,nPerSpf), 4*ones(1,nPerSpf)];

    % If condition has leftover trials, ignore them for SPF analysis
    if size(dHFD,2) > nUsed
        fprintf('Note: using first %d/%d trials for SPF analysis.\n', nUsed, size(dHFD,2));
    end

    % Electrode display list
    NONHP_show = {'Cz','C4','Fz','F3'};
    all_show   = [HP_order, NONHP_show];
    nHP_show   = numel(HP_order);
    nAll_show  = numel(all_show);

    % Colours
    C_BL   = [0.17 0.37 0.54];
    C_ST   = [0.71 0.27 0.11];
    C_HP   = [0.24 0.42 0.32];
    C_GRAY = [0.55 0.55 0.55];

    % ── Figure 1: Main HFD results ───────────────────────────────────────
    fig1 = figure('Name', ['HFD Main Results - ' tag], ...
                  'Position', [50 50 1400 900]);
    fig1.Color = [0.98 0.98 0.96];

    % ── Panel A: Violin BL vs ST ─────────────────────────────────────────
    subplot(2,3,1); hold on;

    for col_i = 1:2
        if col_i == 1
            vals = hp_BL_flat;
            xc = 1;
            col = C_BL;
        else
            vals = hp_ST_flat;
            xc = 2;
            col = C_ST;
        end

        [f_d, xi] = ksdensity(vals, 'NumPoints', 100);
        f_d = f_d / max(f_d) * 0.35;

        fill([xc+f_d, xc-fliplr(f_d)], [xi, fliplr(xi)], col, ...
             'FaceAlpha', 0.6, 'EdgeColor', 'none');

        plot([xc-0.3 xc+0.3], [median(vals) median(vals)], ...
             '-', 'Color', [0 0 0], 'LineWidth', 2.5);
    end

    y_top = max([hp_BL_flat hp_ST_flat]) + 0.005;
    plot([1 1 2 2], [y_top-0.001 y_top y_top y_top-0.001], ...
         'k-', 'LineWidth', 1);

    if p_hp < 0.001
        sig_str = '***';
    elseif p_hp < 0.01
        sig_str = '**';
    elseif p_hp < 0.05
        sig_str = '*';
    else
        sig_str = 'ns';
    end

    text(1.5, y_top+0.001, sprintf('p = %.2e %s', p_hp, sig_str), ...
         'HorizontalAlignment', 'center', 'FontSize', 9);

    set(gca, 'XTick', [1 2], ...
             'XTickLabel', {'Baseline','Stimulus'}, ...
             'FontSize', 10);

    ylabel('HFD', 'FontSize', 10);
    title('A  HP Group: BL vs Stimulus', 'FontSize', 10, 'FontWeight', 'bold');

    text(1, median(hp_BL_flat)-0.003, sprintf('%.4f', median(hp_BL_flat)), ...
         'HorizontalAlignment','center','FontSize',9,'Color',C_BL,'FontWeight','bold');

    text(2, median(hp_ST_flat)-0.003, sprintf('%.4f', median(hp_ST_flat)), ...
         'HorizontalAlignment','center','FontSize',9,'Color',C_ST,'FontWeight','bold');

    box off;

    % ── Panel B: Matched pairs ───────────────────────────────────────────
    subplot(2,3,2); hold on;

    rng(42);
    sample_idx = randperm(numel(hp_BL_flat), min(80, numel(hp_BL_flat)));

    for i = sample_idx
        col = C_ST * 0.8;
        if hp_ST_flat(i) < hp_BL_flat(i)
            col = C_BL * 0.8;
        end

        plot([0 1], [hp_BL_flat(i) hp_ST_flat(i)], ...
             '-', 'Color', [col 0.2], 'LineWidth', 0.7);
    end

    plot([0 1], [median(hp_BL_flat) median(hp_ST_flat)], ...
         'o-', 'Color', [0.75 0.1 0.1], ...
         'LineWidth', 3, 'MarkerSize', 9, ...
         'MarkerFaceColor', [0.75 0.1 0.1]);

    set(gca, 'XTick', [0 1], ...
             'XTickLabel', {'Baseline','Stimulus'}, ...
             'FontSize', 10);

    ylabel('HFD', 'FontSize', 10);

    pct_up = 100 * mean(hp_ST_flat > hp_BL_flat);
    title(sprintf('B  Matched pairs  (%.0f%% ST>BL)', pct_up), ...
          'FontSize', 10, 'FontWeight', 'bold');

    legend('Median', 'Location', 'northwest', 'FontSize', 9);
    box off;

    % ── Panel C: ΔHFD per electrode ──────────────────────────────────────
    subplot(2,3,3); hold on;

    for i = 1:nAll_show
        label = all_show{i};
        ei = find(strcmp(all_labels, label));

        if isempty(ei)
            continue;
        end

        med_d = median(dHFD(ei,:));
        sem_d = std(dHFD(ei,:)) / sqrt(size(dHFD,2));
        pv = pvals(ei);

        col = C_HP;
        if i > nHP_show
            col = C_GRAY;
        end

        bar(i, med_d, 0.6, ...
            'FaceColor', col, ...
            'FaceAlpha', 0.75, ...
            'EdgeColor', [0.1 0.1 0.1], ...
            'LineWidth', 0.5);

        errorbar(i, med_d, sem_d, 'k.', ...
                 'CapSize', 4, 'LineWidth', 1);

        if pv < 0.001
            s = '***';
        elseif pv < 0.01
            s = '**';
        elseif pv < 0.05
            s = '*';
        else
            s = '';
        end

        if ~isempty(s)
            text(i, med_d + 0.001, s, ...
                 'HorizontalAlignment', 'center', ...
                 'FontSize', 8);
        end
    end

    yline(0, 'k--', 'LineWidth', 0.8);
    xline(nHP_show + 0.5, 'k:', 'LineWidth', 0.8, 'Alpha', 0.5);

    set(gca, 'XTick', 1:nAll_show, ...
             'XTickLabel', all_show, ...
             'XTickLabelRotation', 45, ...
             'FontSize', 8);

    ylabel('Median ΔHFD (ST−BL)', 'FontSize', 10);
    title('C  ΔHFD per electrode', 'FontSize', 10, 'FontWeight', 'bold');

    yl = ylim;
    text(nHP_show/2, yl(2)*0.9, 'HP', ...
         'HorizontalAlignment','center','FontSize',8, ...
         'Color',C_HP,'FontAngle','italic');

    text(nHP_show+2.5, yl(2)*0.9, 'Non-HP', ...
         'HorizontalAlignment','center','FontSize',8, ...
         'Color',C_GRAY,'FontAngle','italic');

    box off;

    % ── Panel D: ΔHFD by spatial frequency ───────────────────────────────
    subplot(2,3,4); hold on;

    cols_spf = {C_BL, C_HP, C_ST};

    for si = 1:3
        idx = spf_vals == spf_list(si);

        vals = reshape(dHFD(hp_idx_arr, idx), 1, []);
        bl_v = reshape(HFD_BL(hp_idx_arr, idx), 1, []);
        st_v = reshape(HFD_ST(hp_idx_arr, idx), 1, []);

        pv = signrank(bl_v', st_v');
        m = median(vals);
        se = std(vals) / sqrt(numel(vals));

        bar(si, m, 0.5, ...
            'FaceColor', cols_spf{si}, ...
            'FaceAlpha', 0.75, ...
            'EdgeColor', [0.1 0.1 0.1], ...
            'LineWidth', 0.5);

        errorbar(si, m, se, 'k.', ...
                 'CapSize', 5, 'LineWidth', 1);

        if pv < 0.001
            s = '***';
        elseif pv < 0.01
            s = '**';
        elseif pv < 0.05
            s = '*';
        else
            s = 'ns';
        end

        text(si, m + 0.0006, s, ...
             'HorizontalAlignment', 'center', ...
             'FontSize', 10);
    end

    yline(0, 'k--', 'LineWidth', 0.8);

    set(gca, 'XTick', 1:3, ...
             'XTickLabel', {'1 cpd','2 cpd','4 cpd'}, ...
             'FontSize', 10);

    ylabel('Median ΔHFD (ST−BL)', 'FontSize', 10);
    title('D  ΔHFD by Spatial Frequency', 'FontSize', 10, 'FontWeight', 'bold');

    yl = ylim;
    text(3.45, yl(2)*0.95, sprintf('n=%d/spf', nPerSpf), ...
         'FontSize', 8, ...
         'Color', [0.4 0.4 0.4], ...
         'HorizontalAlignment', 'right');

    box off;

    % ── Panel E: Time-resolved HFD ───────────────────────────────────────
    subplot(2,3,5); hold on;

    fill([tc_time, fliplr(tc_time)], ...
         [tc_med + tc_sem, fliplr(tc_med - tc_sem)], ...
         C_HP, ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none');

    plot(tc_time, tc_med, '-', ...
         'Color', C_HP, ...
         'LineWidth', 2.2);

    xline(0, 'k--', 'LineWidth', 1.2);
    xline(0.8, 'k:', 'LineWidth', 0.8);

    xregion(-0.5, 0, ...
            'FaceColor', C_BL, ...
            'FaceAlpha', 0.06, ...
            'HandleVisibility', 'off');

    xregion(0.25, 0.75, ...
            'FaceColor', C_ST, ...
            'FaceAlpha', 0.07, ...
            'HandleVisibility', 'off');

    xlabel('Time (s)', 'FontSize', 10);
    ylabel('HFD (kmax=5)', 'FontSize', 10);
    title('E  Time-resolved HFD — HP group', 'FontSize', 10, 'FontWeight', 'bold');

    text(0.02, min(tc_med)-0.001, 'Stimulus onset', ...
         'FontSize', 8, ...
         'Rotation', 90);

    box off;

    % ── Panel F: Power spectrum ──────────────────────────────────────────
    subplot(2,3,6); hold on;

    semilogy(f, psd_bl, '-', ...
             'Color', C_BL, ...
             'LineWidth', 1.8, ...
             'DisplayName', 'Baseline');

    semilogy(f, psd_st, '-', ...
             'Color', C_ST, ...
             'LineWidth', 1.8, ...
             'DisplayName', 'Stimulus');

    xregion(8, 12, ...
            'FaceColor', [0.2 0.7 0.2], ...
            'FaceAlpha', 0.12, ...
            'HandleVisibility', 'off');

    xregion(30, 70, ...
            'FaceColor', [0.8 0.2 0.2], ...
            'FaceAlpha', 0.08, ...
            'HandleVisibility', 'off');

    text(10, max(psd_bl)*0.35, '\alpha', ...
         'FontSize', 14, ...
         'Color', [0.1 0.6 0.1], ...
         'HorizontalAlignment', 'center');

    text(50, max(psd_bl)*0.04, '\gamma', ...
         'FontSize', 14, ...
         'Color', [0.7 0.1 0.1], ...
         'HorizontalAlignment', 'center');

    xlabel('Frequency (Hz)', 'FontSize', 10);
    ylabel('PSD (\muV^2/Hz)', 'FontSize', 10);
    title('F  Power Spectrum — Oz electrode', 'FontSize', 10, 'FontWeight', 'bold');

    legend('Location', 'northeast', 'FontSize', 9);
    box off;

    % ── Global title ─────────────────────────────────────────────────────
    sgtitle({'G1 Project — Preliminary Results', ...
         sprintf('Subject %s · %s · HP Electrodes  (n=%d trials, Fs=%d Hz)', ...
                 subjectID, tag, nTrials, Fs)}, ...
         'FontSize', 12, ...
         'FontWeight', 'bold');
    % ── Save condition-specific figures ─────────────────────────────────
    figFile = fullfile(outputDir, ['Fig1_HFD_Main_' subjectID '_' tag '.fig']);
    pngFile = fullfile(outputDir, ['Fig1_HFD_Main_' subjectID '_' tag '.png']);
    saveas(fig1, figFile);
    saveas(fig1, pngFile);

    fprintf('Figure saved:\n  %s\n  %s\n', figFile, pngFile);

    % ── Summary table printed to command window ─────────────────────────
    fprintf('\n%s\n', repmat('=',1,72));
    fprintf('RESULTS TABLE — Subject %s, %s  (Fs=%d Hz, n=%d trials)\n', ...
        subjectID, tag, Fs, nTrials);
    fprintf('%s\n', repmat('=',1,72));
    fprintf('%-8s %7s %7s %8s %12s %5s\n', ...
            'Electrode','BL','ST','ΔHFD','p-value','sig');
    fprintf('%s\n', repmat('-',1,72));

    for i = 1:numel(all_labels)
        label = all_labels{i};

        if ~ismember(label, all_show)
            continue;
        end

        i_in_all = find(strcmp(all_labels, label));
        pv = pvals(i_in_all);

        if pv < 0.001
            s = '***';
        elseif pv < 0.01
            s = '**';
        elseif pv < 0.05
            s = '*';
        else
            s = 'ns';
        end

        hp_tag = '    ';
        if ismember(label, HP_order)
            hp_tag = '[HP]';
        end

        fprintf('%s%-5s %7.4f %7.4f %+8.4f %12.2e %5s\n', ...
                hp_tag, label, ...
                median(HFD_BL(i_in_all,:)), ...
                median(HFD_ST(i_in_all,:)), ...
                median(dHFD(i_in_all,:)), ...
                pv, s);
    end

    fprintf('%s\n', repmat('-',1,72));

    if p_hp < 0.001
        hp_sig = '***';
    elseif p_hp < 0.01
        hp_sig = '**';
    elseif p_hp < 0.05
        hp_sig = '*';
    else
        hp_sig = 'ns';
    end

    fprintf('%-13s %7.4f %7.4f %+8.4f %12.2e %5s\n', ...
            'HP pooled', ...
            median(hp_BL_flat), ...
            median(hp_ST_flat), ...
            median(hp_ST_flat - hp_BL_flat), ...
            p_hp, hp_sig);

    fprintf('Finished plotting: %s\n', tag);

    close(fig1);

end

fprintf('\nBatch plotting complete.\n');