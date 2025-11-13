%% ========================================================================
%  EXPERIMENT 1: MASTER ANALYSIS SCRIPT
%  - Adaptation (Phase 2)          -> counterbalanced group is flipped;
%  change is hand angle wrt error on previous trial
%  - Temporal Context Effects      -> counterbalanced group is flipped;
%  change in hand angle wrt error on the SAME trial 
%  - Generalization (Phase 3)      -> handAng; counterbalanced subgroups
%  remain as is for better visualization
% ========================================================================

clear; clc; close all;
set(0, 'DefaultFigureUnits', 'pixels', ...
       'DefaultFigurePosition', [200 200 600 420]);


%% FILES
data_dir  = fullfile('data');  % or '.' or '../data', etc.
exp_file  = fullfile(data_dir, 'twocontext_data.csv');
ctrl_file = fullfile(data_dir, 'control_data.csv');

% Counterbalanced (flip) participant IDs for EXPERIMENTAL GROUP
flip_ids = [2035 2040 2070 2092 2160 2183 2184 2185 2186];

% Control exclusions for generalization
ctrl_exclude_ids = [2234 2249];

%% ======================= LOAD BASE DATA (UNFLIPPED) =====================

data_raw = readtable(exp_file);


data = data_raw;

%% ======================= APPLY FLIPS (ADAPTATION/TEMPORAL) ==============
% Flip signs for counterbalanced participants in the EXPERIMENTAL data copy
% (This copy will be used ONLY for adaptation & temporal context effects.)

flip_idx = ismember(data.id, flip_ids);

% 1) Rotation directions
if ismember('newRotDir', data.Properties.VariableNames)
    data.newRotDir(flip_idx) = -data.newRotDir(flip_idx);
end
if ismember('rotDir', data.Properties.VariableNames)
    data.rotDir(flip_idx) = -data.rotDir(flip_idx);
end

% 2) Behavioral angles (only if they exist)
varsToFlip = {'handAng', 'headAng', 'velAng', 'cursorAng', 'newAdaptationAng'};
for i = 1:numel(varsToFlip)
    var = varsToFlip{i};
    if ismember(var, data.Properties.VariableNames)
        data.(var)(flip_idx) = -data.(var)(flip_idx);
    end
end

%% ========================================================================
%  PART 1: ADAPTATION BY PREVIOUS ROTATION DIRECTION (PHASE 2, newRotDir, newAdaptationAng)
% ========================================================================

phase2Data = data(data.phase == 2, :);

% Match the text label: use trials 21–420 within Phase 2
trialData = phase2Data(phase2Data.subTrial >= 21 & phase2Data.subTrial <= 420, :);

% Focus only on rotation directions -1 and +1
uniqueRotDirs = [-1 1];

participants = unique(trialData.id);
numParticipants = numel(participants);
numRotDirs = numel(uniqueRotDirs);

% Matrix: participant × rotDir mean newAdaptationAng
avgAdaptByRot = nan(numParticipants, numRotDirs);

for p = 1:numParticipants
    pid = participants(p);
    participantData = trialData(trialData.id == pid, :);

    for r = 1:numRotDirs
        rotDir = uniqueRotDirs(r);
        rotData = participantData(participantData.newRotDir == rotDir, :);

        if ~isempty(rotData)
            avgAdaptByRot(p, r) = mean(rotData.newAdaptationAng, 'omitnan');
        end
    end
end

% Group stats
groupMeans = mean(avgAdaptByRot, 1, 'omitnan');
groupSEM   = std(avgAdaptByRot, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(avgAdaptByRot), 1));
groupSD    = std(avgAdaptByRot, 0, 1, 'omitnan');

% -------------------- Plot: Adaptation by Rotation Direction -------------
figure('Name', 'Adaptation by Rotation Direction');

hold on;
barWidth = 0.6;
bar(uniqueRotDirs, groupMeans, barWidth, 'FaceColor', [0.4 0.6 0.8]);
errorbar(uniqueRotDirs, groupMeans, groupSEM, 'k', 'LineStyle', 'none', ...
    'LineWidth', 1.5, 'CapSize', 10);

% Individual participant dots
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
title('Adaptation Effect', ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on; box off;
ylim([-3 3]);
yline(0, '--k', 'LineWidth', 1);

xticks(uniqueRotDirs);
xticklabels({'-1', '+1'});

% -------------------- Print Summary & T-tests ----------------------------
fprintf('\n=== PART 1: ADAPTATION  ===\n');
fprintf('Delta hand angles by previous error direction (newRotDir):\n');
for r = 1:numRotDirs
    fprintf('RotDir %.1f: %.2f ± %.2f (mean ± SEM), SD = %.2f\n', ...
        uniqueRotDirs(r), groupMeans(r), groupSEM(r), groupSD(r));
end

% 1) rotDir = -1 vs 0
rotNeg1Data = avgAdaptByRot(:, 1);
rotNeg1Data = rotNeg1Data(~isnan(rotNeg1Data));

fprintf('\n--- T-Test Results (against zero) ---\n');
if ~isempty(rotNeg1Data)
    [~, p_neg1, ~, stats_neg1] = ttest(rotNeg1Data, 0);
    fprintf('RotDir = -1: mean = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        mean(rotNeg1Data), stats_neg1.df, stats_neg1.tstat, p_neg1);
else
    fprintf('RotDir = -1: no valid data for t-test.\n');
end

% 2) rotDir = +1 vs 0
rotPos1Data = avgAdaptByRot(:, 2);
rotPos1Data = rotPos1Data(~isnan(rotPos1Data));

if ~isempty(rotPos1Data)
    [~, p_pos1, ~, stats_pos1] = ttest(rotPos1Data, 0);
    fprintf('RotDir = +1: mean = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        mean(rotPos1Data), stats_pos1.df, stats_pos1.tstat, p_pos1);
else
    fprintf('RotDir = +1: no valid data for t-test.\n');
end

% 3) Paired comparison +1 vs -1
validParticipants = ~isnan(avgAdaptByRot(:, 1)) & ~isnan(avgAdaptByRot(:, 2));
rotNeg1DataPaired = avgAdaptByRot(validParticipants, 1);
rotPos1DataPaired = avgAdaptByRot(validParticipants, 2);

if sum(validParticipants) > 1
    [~, p_pair, ~, stats_pair] = ttest(rotPos1DataPaired, rotNeg1DataPaired);
    mean_diff = mean(rotPos1DataPaired - rotNeg1DataPaired);
    fprintf('Paired +1 vs -1: mean diff = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        mean_diff, stats_pair.df, stats_pair.tstat, p_pair);
else
    fprintf('recheck');
end



%% ========================================================================
%  PART 2: TEMPORAL CONTEXT EFFECTS- hand angle change ON SAME TRIAL
%          (Phase 2,rotDir, newAdaptationAng, with flips)
% ========================================================================

fprintf('\n\n=============================================================\n');
fprintf('PART 2: TEMPORAL CONTEXT EFFECTS\n');
fprintf('=============================================================\n');

phase2 = data(data.phase == 2, :);
participants = unique(phase2.id);

% We now use ONLY newAdaptationAng for temporal effects
measures = {'newAdaptationAng'};
phase_names = {'early', 'late'};
trial_ranges = [21 220; 221 420];  % early / late windows

results = struct();

for m = 1:numel(measures)
    var = measures{m};
    for p = 1:2
        range = trial_ranges(p, :);
        phaseData = phase2(phase2.subTrial >= range(1) & phase2.subTrial <= range(2), :);

        for d = [-1 1]
            condData = phaseData(phaseData.rotDir == d, :);
            fieldname = sprintf('rot_%s', strrep(num2str(d), '-', 'm'));

            for i = 1:numel(participants)
                pid = participants(i);
                subjData = condData(condData.id == pid, :);
                results.(var).(phase_names{p}).(fieldname)(i, 1) = ...
                    mean(subjData.(var), 'omitnan');
            end
        end
    end

    % Condition averages
    earlyPlus = results.(var).early.rot_1;
    earlyMinus = results.(var).early.rot_m1;
    latePlus  = results.(var).late.rot_1;
    lateMinus = results.(var).late.rot_m1;

    earlyDiff = earlyMinus - earlyPlus; % (error -1 minus +1)
    lateDiff  = lateMinus - latePlus;
    diffDiff  = lateDiff - earlyDiff;

    results.(var).earlyDiff = earlyDiff;
    results.(var).lateDiff  = lateDiff;
    results.(var).diffDiff  = diffDiff;

    % ------------------- T-tests (participant-level diffs) ---------------
    fprintf('\n=== PART 1: CONTEXT EFFECT  ===\n');
    [~, p1, ci1, stats1] = ttest(earlyDiff);
    [~, p2, ci2, stats2] = ttest(lateDiff);
    [~, p3, ci3, stats3] = ttest(diffDiff);

    fprintf('Early Diff (-1 - +1): t(%d)=%.3f, p=%.4f, CI=[%.2f %.2f]\n', ...
        stats1.df, stats1.tstat, p1, ci1(1), ci1(2));
    fprintf('Late Diff  (-1 - +1): t(%d)=%.3f, p=%.4f, CI=[%.2f %.2f]\n', ...
        stats2.df, stats2.tstat, p2, ci2(1), ci2(2));
    fprintf('Late - Early Diff:     t(%d)=%.3f, p=%.4f, CI=[%.2f %.2f]\n', ...
        stats3.df, stats3.tstat, p3, ci3(1), ci3(2));

    % ------------------- Plot: Early/Late/Diff of Diffs -----------------
    figure; hold on;
    means = [mean(earlyDiff,'omitnan'), ...
             mean(lateDiff,'omitnan'), ...
             mean(diffDiff,'omitnan')];

    n_subj = sum(~isnan(earlyDiff)); % assume same N for all three, roughly
    sems = [ std(earlyDiff,0,'omitnan')/sqrt(n_subj), ...
             std(lateDiff,0,'omitnan')/sqrt(n_subj), ...
             std(diffDiff,0,'omitnan')/sqrt(n_subj) ];

    bar(1:3, means, 'FaceColor', [0.6 0.6 0.6]);
    errorbar(1:3, means, sems, '.k', 'LineWidth', 1.5);
    title('Context Effect', ...
    'FontSize', 14, 'FontWeight', 'bold');

    scatter(ones(size(earlyDiff)), earlyDiff, 40, 'b', 'filled');
    scatter(2*ones(size(lateDiff)), lateDiff, 40, 'r', 'filled');
    scatter(3*ones(size(diffDiff)), diffDiff, 40, 'k', 'filled');

    xticks(1:3);
    xticklabels({'Early', 'Late', 'Late - Early'});
    ylabel(sprintf('Difference (-1 minus +1)\n[%s]', var));
    xlim([0.5 3.5]);
    ylim([-1 4]);
    box off;
end


%% ========================================================================
%  PART 3: GENERALIZATION (PHASE 3, handAng, NO FLIPS)
% ========================================================================

fprintf('\n\n=============================================================\n');
fprintf('PART 3: GENERALIZATION (Phase 3, hand Angle, NO flips)\n');
fprintf('=============================================================\n');

goDelays = [500 750 1000 1250 1500];
meas = 'handAng'; % absolute hand angle

% -------- CONTROL GROUP (Phase 3, exclusions, mean-centered per subject) -
C = readtable(ctrl_file);
C = C(C.phase == 3, :);
C = C(~ismember(C.id, ctrl_exclude_ids), :);

for pid = unique(C.id)'
    idx = C.id == pid;
    C.(meas)(idx) = C.(meas)(idx) - mean(C.(meas)(idx), 'omitnan');
end

% -------- EXPERIMENTAL GROUP (Phase 3, NO flips, mean-centered) ----------
E = data_raw;                 % use UNFLIPPED copy
E = E(E.phase == 3, :);

for pid = unique(E.id)'
    idx = E.id == pid;
    E.(meas)(idx) = E.(meas)(idx) - mean(E.(meas)(idx), 'omitnan');
end

% Split experimental into two context subgroups
% Subgroup 1: non-flipped IDs (e.g., CW/CCW)
% Subgroup 2: flipped IDs   (e.g., CCW/CW)
E_sub1 = E(~ismember(E.id, flip_ids), :);   % Subgroup 1
E_sub2 = E( ismember(E.id, flip_ids), :);   % Subgroup 2



% Helper: subject × delay matrix
compute_group_matrix = @(T) arrayfun(@(pid) ...
    arrayfun(@(dly) mean(T.(meas)(T.id==pid & T.goDelay==dly), 'omitnan'), goDelays), ...
    unique(T.id)', 'UniformOutput', false);
group_mean_matrix = @(T) cell2mat(compute_group_matrix(T)');  % numeric

ctrl_mat   = group_mean_matrix(C);
sub1_mat   = group_mean_matrix(E_sub1);
sub2_mat   = group_mean_matrix(E_sub2);

% Group means & SEs
ctrl_mean = mean(ctrl_mat, 1, 'omitnan');
ctrl_se   = std(ctrl_mat, 0, 1, 'omitnan') ./ sqrt(size(ctrl_mat, 1));

sub1_mean = mean(sub1_mat, 1, 'omitnan');
sub1_se   = std(sub1_mat, 0, 1, 'omitnan') ./ sqrt(size(sub1_mat, 1));

sub2_mean = mean(sub2_mat, 1, 'omitnan');
sub2_se   = std(sub2_mat, 0, 1, 'omitnan') ./ sqrt(size(sub2_mat, 1));

% ------------------------ Plot: Generalization ---------------------------
figure('Color','w','Name','Phase 3: handAng by Delay (no flips)');
hold on; grid on; box on;

errorbar(goDelays, sub1_mean, sub1_se, '-s', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Context group - Subgroup 1 (CW/CCW)');
errorbar(goDelays, sub2_mean, sub2_se, '-^', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'Context group - Subgroup 2 (CCW/CW)');

xlabel('Planning Delay (ms)', 'FontSize', 12);
ylabel('Mean-centered handAng (deg)', 'FontSize', 12);
legend('Location','best', 'FontSize', 11);
title('Phase 3 Generalization: Mean-centered hand Ang vs Delay', 'FontSize', 13);
set(gca, 'FontSize', 11);

% ------------------------ Mixed ANOVAs -----------------------------------

% Build wide tables: columns d500,...,d1500, plus id, group
delays = goDelays;

ctrl_tbl = array2table(ctrl_mat, 'VariableNames', compose('d%d', delays));
ctrl_tbl.id    = unique(C.id);
ctrl_tbl.group = repmat({'Control'}, height(ctrl_tbl), 1);

sub1_tbl = array2table(sub1_mat, 'VariableNames', compose('d%d', delays));
sub1_tbl.id    = unique(E_sub1.id);
sub1_tbl.group = repmat({'Subgroup1'}, height(sub1_tbl), 1);

sub2_tbl = array2table(sub2_mat, 'VariableNames', compose('d%d', delays));
sub2_tbl.id    = unique(E_sub2.id);
sub2_tbl.group = repmat({'Subgroup2'}, height(sub2_tbl), 1);

within = table(categorical(goDelays'), 'VariableNames', {'Delay'});

% Control vs Subgroup 1
tbl1 = [ctrl_tbl; sub1_tbl];
rm1  = fitrm(tbl1, 'd500-d1500 ~ group', 'WithinDesign', within);
anova1 = ranova(rm1, 'WithinModel', 'Delay');

fprintf('\nControl vs Subgroup 1 (CW/CCW)\n');
disp(anova1);

% Control vs Subgroup 2
tbl2 = [ctrl_tbl; sub2_tbl];
rm2  = fitrm(tbl2, 'd500-d1500 ~ group', 'WithinDesign', within);
anova2 = ranova(rm2, 'WithinModel', 'Delay');

fprintf('\nControl vs Subgroup 2 (CCW/CW)\n');
disp(anova2);

fprintf('COMPLETE\n');
