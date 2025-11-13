%  THREE-DELAY EXPERIMENT: ANALYSIS SCRIPT
%  GRAPHS:
%   1) Phase 2 adaptation:- change (delta) in hand angle by previous error direction (newRotDir = -15 vs +15)
%
%   2) Phase 2 temporal effect 
%        - Early vs late (500 − 1000)
%        - Early vs late (1500 − 1000)
%
%   3) Phase 3 generalization:
%        - handAng, mean-centered per participant, across
%          goDelay = [500 750 1000 1250 1500]
%
%   4) Model comparison on generalization (using 5-delay averages):
%        - Quadratic vs linear


clear; clc; close all;

set(0, 'DefaultFigureUnits',   'pixels', ...
       'DefaultFigurePosition',[200 200 650 450]);

data_dir  = fullfile('data');
exp_file  = fullfile(data_dir, 'threecontext_data.csv');

data_raw = readtable(exp_file);

% Participant groups (counterbalanced)
group1_ids = [2300 2278 2279 2280 2281 2282 2283 2284 2285 2286];
group2_ids = [2304 2305 2308 2325 2326 2328 2329 2331 2332 2327];  


% Flip signs for group 2 so that clamp/error directions are aligned
% across subjects. This flipped table is used for ALL analyses.

data   = data_raw;
idx_g2 = ismember(data.id, group2_ids);

% Rotation directions
if ismember('rotDir', data.Properties.VariableNames)
    data.rotDir(idx_g2) = -data.rotDir(idx_g2);
end
if ismember('newRotDir', data.Properties.VariableNames)
    data.newRotDir(idx_g2) = -data.newRotDir(idx_g2);
end

% Behavioral angles
flipVars = {'handAng','headAng','velAng','cursorAng','newAdaptationAng'};
for i = 1:numel(flipVars)
    v = flipVars{i};
    if ismember(v, data.Properties.VariableNames)
        data.(v)(idx_g2) = -data.(v)(idx_g2);
    end
end

fprintf('Applied sign flip to group 2 for rotDir/newRotDir and angle variables.\n');

%  PART 1: PHASE 2 ADAPTATION BY PREVIOUS ERROR (newRotDir - newAdaptationAng)

fprintf('PART 1: PHASE 2 ADAPTATION BY PREVIOUS ERROR\n');

valid_ids = [group1_ids group2_ids];
phase2 = data(data.phase == 2 & ismember(data.id, valid_ids), :);

phase2 = phase2(phase2.subTrial >= 21 & phase2.subTrial <= 420, :);

participants = unique(phase2.id);
nSubj        = numel(participants);
uniqueRotDirs = [-1 1];

avgAdaptByRot = nan(nSubj, 2);

for i = 1:nSubj
    pid  = participants(i);
    pdat = phase2(phase2.id == pid, :);
    for r = 1:2
        d    = uniqueRotDirs(r);
        tr   = pdat(pdat.newRotDir == d, :);
        avgAdaptByRot(i,r) = mean(tr.newAdaptationAng, 'omitnan');
    end
end

% group summary
groupMean = mean(avgAdaptByRot, 1, 'omitnan');
groupSE   = std(avgAdaptByRot, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(avgAdaptByRot),1));

% ---------- Plot 
figure('Name','Phase 2: Adaptation by Previous Error (newRotDir)');
hold on;

barWidth = 0.4;
bar(uniqueRotDirs, groupMean, barWidth, 'FaceColor', [0.4 0.6 0.8]);
errorbar(uniqueRotDirs, groupMean, groupSE, 'k', ...
    'LineStyle','none', 'LineWidth',1.5, 'CapSize',10);

xlabel('Error direction on trial t-1 (newRotDir)');
ylabel('Δ hand angle (deg) [newAdaptationAng]');
title('Three-delay: Phase 2 Adaptation by Previous Error');
yline(0,'--k','LineWidth',1);

yMax = max(abs(groupMean) + groupSE);
if yMax < 1e-6, yMax = 1; end
yPad = 0.15 * yMax;
ylim([-yMax - yPad, yMax + yPad]);

xlim([-1.6 1.6]);
set(gca,'XTick',uniqueRotDirs,'XTickLabel',{'-1','+1'}, ...
    'LineWidth',1.5,'FontSize',12,'TickDir','out','Box','off');
grid on;

% ---------- T-tests 
rotNeg = avgAdaptByRot(:,1); rotNeg = rotNeg(isfinite(rotNeg));
rotPos = avgAdaptByRot(:,2); rotPos = rotPos(isfinite(rotPos));

fprintf(' Adaptation statistics ');

if ~isempty(rotNeg)
    [~, p_neg, ~, s_neg] = ttest(rotNeg, 0);
    fprintf('newRotDir = -1: mean = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        mean(rotNeg), s_neg.df, s_neg.tstat, p_neg);
end

if ~isempty(rotPos)
    [~, p_pos, ~, s_pos] = ttest(rotPos, 0);
    fprintf('newRotDir = +1: mean = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        mean(rotPos), s_pos.df, s_pos.tstat, p_pos);
end

valid = isfinite(avgAdaptByRot(:,1)) & isfinite(avgAdaptByRot(:,2));
if sum(valid) > 1
    [~, p_pair, ~, s_pair] = ttest(avgAdaptByRot(valid,2), avgAdaptByRot(valid,1));
    diffMean = mean(avgAdaptByRot(valid,2) - avgAdaptByRot(valid,1));
    fprintf('+1 vs -1 (paired): mean diff = %.3f, t(%d) = %.3f, p = %.4f\n', ...
        diffMean, s_pair.df, s_pair.tstat, p_pair);
else
    fprintf('recheck');
end


%  PART 2: PHASE 2 TEMPORAL EFFECTS ( 500–1000 & 1500–1000)

fprintf('TEMPORAL EFFECTS');

% early / late windows
earlyRange = [21 220];
lateRange  = [221 420];
goDelays   = [500 1000 1500];

participants2 = unique(phase2.id);
nSubj2        = numel(participants2);

% earlyLateMeans
earlyLateMeans = nan(nSubj2, 6);

for i = 1:nSubj2
    pid   = participants2(i);
    pdata = phase2(phase2.id == pid, :);

    for d = 1:numel(goDelays)
        delay = goDelays(d);

        early = pdata(pdata.goDelay == delay & ...
                      pdata.subTrial >= earlyRange(1) & pdata.subTrial <= earlyRange(2), :);
        late  = pdata(pdata.goDelay == delay & ...
                      pdata.subTrial >= lateRange(1)  & pdata.subTrial <= lateRange(2), :);

        earlyLateMeans(i, (d-1)*2 + 1) = mean(early.newAdaptationAng, 'omitnan');
        earlyLateMeans(i, (d-1)*2 + 2) = mean(late.newAdaptationAng,  'omitnan');
    end
end

%% make temporal-effect figure
makeTemporalFig = @(earlyDiff, lateDiff, nameTag) ...
    localTemporalFigure(earlyDiff, lateDiff, nameTag);

% --------------- 500 − 1000 comparison 

early_500  = earlyLateMeans(:,1);
late_500   = earlyLateMeans(:,2);
early_1000 = earlyLateMeans(:,3);
late_1000  = earlyLateMeans(:,4);

% (500 − 1000) for early and late
earlyDiff_500_1000 = early_500 - early_1000;
lateDiff_500_1000  = late_500  - late_1000;

makeTemporalFig(earlyDiff_500_1000, lateDiff_500_1000, 'Δ(500−1000)');

% --------------- 1500 − 1000 comparison 
early_1500 = earlyLateMeans(:,5);
late_1500  = earlyLateMeans(:,6);

% (1500 − 1000) for early and late
earlyDiff_1500_1000 = early_1500 - early_1000;
lateDiff_1500_1000  = late_1500  - late_1000;

makeTemporalFig(earlyDiff_1500_1000, lateDiff_1500_1000, 'Δ(1500−1000)');

% PHASE 3 GENERALIZATION (handAng, mean-centered)

fprintf('GENERALIZATION ');

goDelaysP3 = [500 750 1000 1250 1500];
meas       = 'handAng';

valid_ids = [group1_ids group2_ids];
phase3 = data(data.phase == 3 & ismember(data.id, valid_ids), :);
participants3 = unique(phase3.id);
nSubj3        = numel(participants3);


% mean-center handAng per participant (trial level)
for i = 1:nSubj3
    pid = participants3(i);
    idx = phase3.id == pid;
    phase3.(meas)(idx) = phase3.(meas)(idx) - mean(phase3.(meas)(idx), 'omitnan');
end

% participant × delay means (5-value per participant)
nDelay   = numel(goDelaysP3);
subjDelayMeans = nan(nSubj3, nDelay);

for i = 1:nSubj3
    pid   = participants3(i);
    pdata = phase3(phase3.id == pid, :);
    for d = 1:nDelay
        delay = goDelaysP3(d);
        dd    = pdata(pdata.goDelay == delay, :);
        subjDelayMeans(i,d) = mean(dd.(meas), 'omitnan');
    end
end

groupMeanP3 = mean(subjDelayMeans, 1, 'omitnan');
groupSEP3   = std(subjDelayMeans, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(subjDelayMeans),1));

% ---------- Generalization plot 
figure('Name','Phase 3: Generalization (handAng)');
hold on;

errorbar(goDelaysP3, groupMeanP3, groupSEP3, '-o', ...
    'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.2 0.4 0.6], 'Color', [0.2 0.4 0.6], ...
    'CapSize', 10);

xlabel('Go Delay (ms)');
ylabel('Mean-centered handAng (deg)');
title('Phase 3: Generalization across planning delay');

yAll = [groupMeanP3 + groupSEP3, groupMeanP3 - groupSEP3];
yMax = max(abs(yAll));
if yMax < 1e-6, yMax = 1; end
yPad = 0.15 * yMax;
ylim([-yMax - yPad, yMax + yPad]);

xMargin = 50;
xlim([min(goDelaysP3)-xMargin, max(goDelaysP3)+xMargin]);

set(gca,'LineWidth',1.5,'FontSize',12,'TickDir','out','Box','on');
grid on;

%   MODEL COMPARISON 

fprintf('PART 4: MODEL COMPARISON ON GENERALIZATION (QUAD vs LIN)\n');

% Use the 5 averaged values per participant (subjDelayMeans)
gen = subjDelayMeans;              % nSubj3 x 5
[nsubs, nprobes] = size(gen);

% Normalize (mean-center per subject based on these 5 values)
gen = gen - mean(gen, 2, 'omitnan');

% ----- Preallocate -----
r2       = nan(nsubs,1);
sse      = nan(nsubs,1);
aic      = nan(nsubs,1);
r2_lin   = nan(nsubs,1);
sse_lin  = nan(nsubs,1);
aic_lin  = nan(nsubs,1);
param    = nan(nsubs,3);  

x = (1:nprobes)';   % index space [1..5]
np = numel(x);

for j = 1:nsubs
    y = gen(j,:)';

    % Quadratic fit (poly2) with nonnegative curvature (a >= 0)
    fo = fitoptions('poly2', 'Lower', [0, -Inf, -Inf]);
    try
        [f, gof] = fit(x, y, 'poly2', fo);
        r2(j)  = gof.rsquare;
        sse(j) = gof.sse;
        kp     = 3;
        aic(j) = np * log(sse(j)/np) + 2*kp;

        param(j,1) = f.p1;  % a
        param(j,2) = f.p2;  % b
        param(j,3) = f.p3;  % c
    catch
        warning('Subject %d: quadratic fit failed.', j);
    end

    % Linear fit (poly1)
    try
        [f_lin, gof_lin] = fit(x, y, 'poly1');
        r2_lin(j)   = gof_lin.rsquare;
        sse_lin(j)  = gof_lin.sse;
        kp_lin      = 2;
        aic_lin(j)  = np * log(sse_lin(j)/np) + 2*kp_lin;
    catch
        warning('Subject %d: linear fit failed.', j);
    end
end

% Basic AIC 
n      = nprobes;   % 5 probe delays
k_lin  = 2;         % linear params
k_quad = 3;         % quadratic params;

AIC_lin  = n*log(sse_lin./n) + 2*k_lin;
AIC_quad = n*log(sse./n)     + 2*k_quad;

%% ----- Plot: Rsq and AIC comparison 
figure('Name','Model comparison: R^2 and AIC');
tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

% (1) Rsq comparison
nexttile; hold on; box on;
meanR2_quad = mean(r2,    'omitnan');
meanR2_lin  = mean(r2_lin,'omitnan');
seR2_quad   = std(r2,    'omitnan') / sqrt(sum(isfinite(r2)));
seR2_lin    = std(r2_lin,'omitnan') / sqrt(sum(isfinite(r2_lin)));

barWidth = 0.45;
bar(1, meanR2_quad, barWidth);
bar(2, meanR2_lin,  barWidth);
errorbar([1 2], [meanR2_quad meanR2_lin], [seR2_quad seR2_lin], ...
         'k', 'LineStyle','none', 'LineWidth',1.5, 'CapSize',10);

set(gca,'XTick',[1 2],'XTickLabel',{'Quadratic','Linear'}, ...
    'LineWidth',1.5,'FontSize',12,'TickDir','out');
ylabel('R^2');
title('Model R^2 (mean \pm SE)');

yAllR2 = [meanR2_quad+seR2_quad, meanR2_lin+seR2_lin];
ylim([0, max(yAllR2)*1.2]);

% (2) AIC sum comparison
nexttile; hold on; box on;

x12 = [1 2];

% per-subject AIC in index space 
sumAIC_quad_points = nprobes*log(sse./nprobes) + 2*3;
sumAIC_lin_points  = nprobes*log(sse_lin./nprobes) + 2*2;

% plot the summed AIC across subjects 
plot(x12, [nansum(sumAIC_quad_points) nansum(sumAIC_lin_points)], 'o-', ...
     'LineWidth', 2, 'MarkerFaceColor', 'k');

set(gca,'XTick',x12,'XTickLabel',{'Quadratic','Linear'}, ...
    'LineWidth',2,'FontSize',11);
ylabel('AIC');
title('Quadratic vs. Linear (points)');
grid on;


%% ----- Quadratic fit on group-mean generalization curve
x_ms  = goDelaysP3(:);                   % [500 750 1000 1250 1500]^T
ybar  = mean(gen, 1, 'omitnan')';        

[f_grp_ms, gof_grp_ms] = fit(x_ms, ybar, 'poly2');   % unconstrained

xMin = min(x_ms);
xMax = max(x_ms);
x_fine = linspace(xMin - 150, xMax + 150, 300)';   
y_fine = feval(f_grp_ms, x_fine);

figure('Name','Quadratic fit to generalization curve');
hold on; box on; grid on;

% Group mean points
plot(x_ms, ybar, 'o', 'LineWidth', 1.5, ...
    'MarkerFaceColor',[0.2 0.4 0.6], 'Color',[0.2 0.4 0.6], ...
    'MarkerSize',8);

% Smoothed quadratic fit
plot(x_fine, y_fine, '-', 'LineWidth', 2, 'Color',[0 0 0]);

xlabel('Probe delay (ms)');
ylabel('Mean-centered HA (deg)');


title(sprintf('Quadratic fit: y = %.3e x^2 + %.3e x + %.3e, R^2 = %.3f', ...
    f_grp_ms.p1, f_grp_ms.p2, f_grp_ms.p3, gof_grp_ms.rsquare));

set(gca,'LineWidth',1.5,'FontSize',12,'TickDir','out');

fprintf('\nQuadratic group-level fit (ms domain):\n');
fprintf('  y = %.6e * x^2 + %.6e * x + %.6e\n', ...
    f_grp_ms.p1, f_grp_ms.p2, f_grp_ms.p3);
fprintf('  R^2 = %.4f\n', gof_grp_ms.rsquare);

fprintf('COMPLETE');

%  temporal-effect plot 

function localTemporalFigure(earlyDiff, lateDiff, nameTag)
 
    validEarly = isfinite(earlyDiff);
    validLate  = isfinite(lateDiff);
    validBoth  = validEarly & validLate;

    earlyDiff = earlyDiff(:);
    lateDiff  = lateDiff(:);

    meanEarly = mean(earlyDiff(validEarly), 'omitnan');
    meanLate  = mean(lateDiff(validLate),   'omitnan');
    meanTE    = mean(lateDiff(validBoth) - earlyDiff(validBoth), 'omitnan');

    seEarly = std(earlyDiff(validEarly), 0, 'omitnan') ./ sqrt(sum(validEarly));
    seLate  = std(lateDiff(validLate),   0, 'omitnan') ./ sqrt(sum(validLate));
    seTE    = std(lateDiff(validBoth) - earlyDiff(validBoth), 0, 'omitnan') ...
              ./ sqrt(sum(validBoth));

    means = [meanEarly, meanLate, meanTE];
    ses   = [seEarly,   seLate,   seTE];

    figure('Name', ['Temporal effect ' nameTag]);
    hold on;

    labels = {['Early '  nameTag], ...
              ['Late '   nameTag], ...
              ['Late−Early ' nameTag]};

    barWidth = 0.4;
    for i = 1:3
        bar(i, means(i), barWidth, 'FaceColor', [0.7 0.7 0.7]);
        errorbar(i, means(i), ses(i), 'k.', 'LineWidth', 1.5, 'CapSize',10);
    end

    % y-limits based on data
    yAll = [means + ses, means - ses];
    yMax = max(abs(yAll));
    if yMax < 1e-6
        yMax = 1; % safeguard
    end
    yPad = 0.15 * yMax;
    ylim([-yMax - yPad, yMax + yPad]);

    xlim([0.5 3.5]);
    ylabel('Δ hand angle (deg)');
    xticks(1:3);
    xticklabels(labels);
    title('Temporal effect');
    yline(0,'--k','LineWidth',1);
    set(gca,'LineWidth',1.5,'FontSize',12,'TickDir','out','Box','off');
    grid on;

    % ----- t-tests -----
    fprintf('\nTemporal effect %s:\n', nameTag);

    if sum(validEarly) > 1
        [~, pE, ~, sE] = ttest(earlyDiff(validEarly), 0);
        fprintf('  Early %s vs 0:      mean = %.3f, t(%d) = %.3f, p = %.4f\n', ...
            nameTag, meanEarly, sE.df, sE.tstat, pE);
    end

    if sum(validLate) > 1
        [~, pL, ~, sL] = ttest(lateDiff(validLate), 0);
        fprintf('  Late %s vs 0:       mean = %.3f, t(%d) = %.3f, p = %.4f\n', ...
            nameTag, meanLate, sL.df, sL.tstat, pL);
    end

    if sum(validBoth) > 1
        [~, pEL, ~, sEL] = ttest(earlyDiff(validBoth), lateDiff(validBoth));
        [~, pTE, ~, sTE] = ttest(lateDiff(validBoth) - earlyDiff(validBoth), 0);
        fprintf('  Late vs Early %s:   t(%d) = %.3f, p = %.4f\n', ...
            nameTag, sEL.df, sEL.tstat, pEL);
        fprintf('  (Late−Early) %s vs 0: t(%d) = %.3f, p = %.4f\n', ...
            nameTag, sTE.df, sTE.tstat, pTE);
    end
end
