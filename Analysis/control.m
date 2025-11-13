%% ========================================================================
%  CONTROL GROUP: ANALYSIS SCRIPT
%
%  - Adaptation (Phase 2)          
%       * no flips
%       * change in hand angle (newAdaptationAng) wrt error direction on
%         previous trial (newRotDir)
%
%  - Temporal Context Effects (Phase 2)
%       * no flips
%       * change in hand angle (newAdaptationAng) wrt error on SAME trial
%         (rotDir), early vs late in Phase 2
%
%  - Generalization (Phase 3)
%       * absolute handAng
%       * repeated-measures ANOVA 


clear; clc; close all;
data_dir  = fullfile('data');  
exp_file  = fullfile(data_dir, 'twocontext_data.csv');
ctrl_file = fullfile(data_dir, 'control_data.csv');


% Make all figures a reasonable uniform size
set(0, 'DefaultFigureUnits',   'pixels', ...
       'DefaultFigurePosition',[200 200 600 420]);

%% ======================= SETTINGS & FILES ===============================

ctrl_exclude_ids = [2234 2249];
goDelays         = [500 750 1000 1250 1500];   % for generalization

% Load full control dataset once
ctrl_raw = readtable(ctrl_file);

%% ========================================================================
%  PART 1: ADAPTATION (PHASE 2, newRotDir, newAdaptationAng)
%          previous-trial error effect
% ========================================================================

fprintf('\n=============================================================\n');
fprintf('PART 1: CONTROL ADAPTATION (Phase 2)\n');
fprintf('=============================================================\n');

% Phase 2 data, apply trial window and exclusions
phase2_all = ctrl_raw(ctrl_raw.phase == 2, :);
phase2     = phase2_all(phase2_all.subTrial >= 21 & phase2_all.subTrial <= 420, :);
phase2     = phase2(~ismember(phase2.id, ctrl_exclude_ids), :);

% Only consider newRotDir = -1 and +1
uniqueRotDirs = [-1 1];

participants   = unique(phase2.id);
numParticipants = numel(participants);
numRotDirs      = numel(uniqueRotDirs);

% Matrix: participant × rotDir mean newAdaptationAng
avgAdaptByRot = nan(numParticipants, numRotDirs);

for p = 1:numParticipants
    pid = participants(p);
    pdata = phase2(phase2.id == pid, :);

    for r = 1:numRotDirs
        rotDir = uniqueRotDirs(r);
        rdata  = pdata(pdata.newRotDir == rotDir, :);
        if ~isempty(rdata)
            avgAdaptByRot(p, r) = mean(rdata.newAdaptationAng, 'omitnan');
        end
    end
end

% Group stats
groupMeans = mean(avgAdaptByRot, 1, 'omitnan');
groupSEM   = std(avgAdaptByRot, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(avgAdaptByRot), 1));
groupSD    = std(avgAdaptByRot, 0, 1, 'omitnan');

% -------------------- Plot: Adaptation by previous error -----------------
figure('Name', 'Control: Adaptation by Previous Error Direction');
hold on;

barWidth = 0.6;
bar(uniqueRotDirs, groupMeans, barWidth, 'FaceColor', [0.4 0.6 0.8]);
errorbar(uniqueRotDirs, groupMeans, groupSEM, 'k', 'LineStyle', 'none', ...
    'LineWidth', 1.5, 'CapSize', 10);

% Individual dots
jitterWidth = 0.1;
for r = 1:numRotDirs
    rotDir = uniqueRotDirs(r);
    for p = 1:numParticipants
        if ~isnan(avgAdaptByRot(p, r))
            scatter(rotDir + (rand-0.5)*jitterWidth, avgAdaptByRot(p, r), ...
                50, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
        end
    end
end

xlabel('Error direction on trial t-1', 'FontSize', 12);
ylabel('Change in HA (deg)', 'FontSize', 12);
title('Control: Adaptation Effect', ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on; box off;
ylim([-2 2]);
yline(0, '--k', 'LineWidth', 1);

xticks(uniqueRotDirs);
xticklabels({'cw', 'ccw'});

% -------------------- Print Summary & T-tests ----------------------------
fprintf('\n=== PART 1: ADAPTATION (CONTROL) ===\n');
fprintf('Delta hand angles by previous error direction:\n');
for r = 1:numRotDirs
    fprintf('RotDir %.1f: %.2f ± %.2f (mean ± SEM), SD = %.2f\n', ...
        uniqueRotDirs(r), groupMeans(r), groupSEM(r), groupSD(r));
end

% T-tests vs zero and between conditions
rotNeg1Data = avgAdaptByRot(:, 1);
rotNeg1Data = rotNeg1Data(~isnan(rotNeg1Data));

rotPos1Data = avgAdaptByRot(:, 2);
rotPos1Data = rotPos1Data(~isnan(rotPos1Data));

fprintf('\n--- T-Test Results (against zero) ---\n');

if ~isempty(rotNeg1Data)
    [~, p_neg1, ~, stats_neg1] = ttest(rotNeg1Data, 0);
    fprintf('RotDir = -1: mean = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        mean(rotNeg1Data), stats_neg1.df, stats_neg1.tstat, p_neg1);
else
    fprintf('recheck');
end

if ~isempty(rotPos1Data)
    [~, p_pos1, ~, stats_pos1] = ttest(rotPos1Data, 0);
    fprintf('RotDir = +1: mean = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        mean(rotPos1Data), stats_pos1.df, stats_pos1.tstat, p_pos1);
else
    fprintf('recheck2');
end

% Paired comparison +1 vs -1
validParticipants = ~isnan(avgAdaptByRot(:, 1)) & ~isnan(avgAdaptByRot(:, 2));
if sum(validParticipants) > 1
    rotNeg1DataPaired = avgAdaptByRot(validParticipants, 1);
    rotPos1DataPaired = avgAdaptByRot(validParticipants, 2);
    [~, p_pair, ~, stats_pair] = ttest(rotPos1DataPaired, rotNeg1DataPaired);
    mean_diff = mean(rotPos1DataPaired - rotNeg1DataPaired);
    fprintf('Paired cw vs ccw: mean diff = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        mean_diff, stats_pair.df, stats_pair.tstat, p_pair);
else
    fprintf('recheck');
end



%% ========================================================================
%  PART 2: TEMPORAL CONTEXT EFFECTS (Phase 2, SAME-TRIAL error)
%          newAdaptationAng, early vs late
% ========================================================================

fprintf('\n\n=============================================================\n');
fprintf('PART 2: CONTROL TEMPORAL CONTEXT EFFECTS');
fprintf('=============================================================\n');

phase2 = phase2_all(~ismember(phase2_all.id, ctrl_exclude_ids), :);
participants = unique(phase2.id);

measures     = {'newAdaptationAng'};   % only what we care about
phase_names  = {'early', 'late'};
trial_ranges = [21 220; 221 420];      % early / late

results = struct();

for m = 1:numel(measures)
    var = measures{m};

    for p = 1:2
        range    = trial_ranges(p, :);
        phaseDat = phase2(phase2.subTrial >= range(1) & phase2.subTrial <= range(2), :);

        for d = [-1 1]
            condDat  = phaseDat(phaseDat.rotDir == d, :);
            fieldname = sprintf('rot_%s', strrep(num2str(d), '-', 'm'));

            for i = 1:numel(participants)
                pid      = participants(i);
                subjData = condDat(condDat.id == pid, :);
                results.(var).(phase_names{p}).(fieldname)(i, 1) = ...
                    mean(subjData.(var), 'omitnan');
            end
        end
    end

    % Condition averages
    earlyPlus  = results.(var).early.rot_1;
    earlyMinus = results.(var).early.rot_m1;
    latePlus   = results.(var).late.rot_1;
    lateMinus  = results.(var).late.rot_m1;

    earlyDiff = earlyMinus - earlyPlus;   % (-1 minus +1)
    lateDiff  = lateMinus  - latePlus;
    diffDiff  = lateDiff   - earlyDiff;

    results.(var).earlyDiff = earlyDiff;
    results.(var).lateDiff  = lateDiff;
    results.(var).diffDiff  = diffDiff;

    % ------------------- T-tests (participant-level diffs) ---------------
    fprintf('\n==== %s (Temporal Context, Control) ====\n', var);
    [~, p1, ci1, stats1] = ttest(earlyDiff);
    [~, p2, ci2, stats2] = ttest(lateDiff);
    [~, p3, ci3, stats3] = ttest(diffDiff);

    fprintf('Early Diff (-1 - +1): t(%d)=%.3f, p=%.4f, CI=[%.2f %.2f]\n', ...
        stats1.df, stats1.tstat, p1, ci1(1), ci1(2));
    fprintf('Late  Diff (-1 - +1): t(%d)=%.3f, p=%.4f, CI=[%.2f %.2f]\n', ...
        stats2.df, stats2.tstat, p2, ci2(1), ci2(2));
    fprintf('Late - Early Diff:    t(%d)=%.3f, p=%.4f, CI=[%.2f %.2f]\n', ...
        stats3.df, stats3.tstat, p3, ci3(1), ci3(2));

    % ------------------- Plot: Early/Late/Diff of Diffs -----------------
    figure; hold on;

    means = [mean(earlyDiff,'omitnan'), ...
             mean(lateDiff,'omitnan'), ...
             mean(diffDiff,'omitnan')];

    n_subj = sum(~isnan(earlyDiff)); % subjects included
    sems   = [std(earlyDiff,0,'omitnan')/sqrt(n_subj), ...
              std(lateDiff,0,'omitnan')/sqrt(n_subj), ...
              std(diffDiff,0,'omitnan')/sqrt(n_subj)];

    bar(1:3, means, 'FaceColor', [0.6 0.6 0.6]);
    errorbar(1:3, means, sems, '.k', 'LineWidth', 1.5);

    scatter(ones(size(earlyDiff)), earlyDiff, 40, 'b', 'filled');
    scatter(2*ones(size(lateDiff)),  lateDiff, 40, 'r', 'filled');
    scatter(3*ones(size(diffDiff)),  diffDiff, 40, 'k', 'filled');

    xticks(1:3);
    xticklabels({'Early', 'Late', 'Late - Early'});
    ylabel('Delta HA (deg) - temporal effect');
    xlabel('Training Phase')
    title('Control: Context Effect');
    xlim([0.5 3.5]);
    ylim([-1 2]);
    box off;
end


%% ========================================================================
%  PART 3: GENERALIZATION (PHASE 3, handAng)
% ========================================================================

fprintf('\n\n=============================================================\n');
fprintf('PART 3: CONTROL GENERALIZATION');
fprintf('=============================================================\n');

meas = 'handAng';

% Phase 3 control data, exclude IDs
phase3 = ctrl_raw(ctrl_raw.phase == 3, :);
phase3 = phase3(~ismember(phase3.id, ctrl_exclude_ids), :);

% Mean-center handAng per participant
for pid = unique(phase3.id)'
    idx = phase3.id == pid;
    phase3.(meas)(idx) = phase3.(meas)(idx) - mean(phase3.(meas)(idx), 'omitnan');
end

% Helper: subject × delay matrix
compute_group_matrix = @(T) arrayfun(@(pid) ...
    arrayfun(@(dly) mean(T.(meas)(T.id==pid & T.goDelay==dly), 'omitnan'), goDelays), ...
    unique(T.id)', 'UniformOutput', false);
group_mean_matrix = @(T) cell2mat(compute_group_matrix(T)');  % numeric

ctrl_mat = group_mean_matrix(phase3);

% Group means & SEs
ctrl_mean = mean(ctrl_mat, 1, 'omitnan');
ctrl_se   = std(ctrl_mat, 0, 1, 'omitnan') ./ sqrt(size(ctrl_mat, 1));

% ------------------------ Plot: Generalization ---------------------------
figure('Color','w','Name','Control Phase 3: handAng by Delay');
hold on; grid on; box on;

errorbar(goDelays, ctrl_mean, ctrl_se, '-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Control');

xlabel('Planning Delay (ms)', 'FontSize', 12);
ylabel('Mean-centered handAng (deg)', 'FontSize', 12);
ylim([-1.5 1.5])
legend('Location','best', 'FontSize', 11);
title('Control Group: Phase 3 Generalization (handAng vs Delay)', 'FontSize', 13);
set(gca, 'FontSize', 11);

fprintf('\nControl participants (Phase 3, after exclusions): %d\n', numel(unique(phase3.id)));

% ------------------------ Repeated-Measures ANOVA -----------------------

delays = goDelays;
ctrl_tbl = array2table(ctrl_mat, 'VariableNames', compose('d%d', delays));
ctrl_tbl.id = unique(phase3.id);

within = table(categorical(goDelays'), 'VariableNames', {'Delay'});

% model with only within-subjects factor Delay
rm_ctrl = fitrm(ctrl_tbl, 'd500-d1500 ~ 1', 'WithinDesign', within);
anova_ctrl = ranova(rm_ctrl, 'WithinModel', 'Delay');

fprintf('\n=== Control Group RM ANOVA (handAng across Delay) ===\n');
disp(anova_ctrl);

fprintf(' COMPLETE');
