function run_fpp_ode_task()
% Run one fixed-point-preserving nonlinear ODE trajectory.

taskId = local_env_number('FPP_TASK_ID', str2double(getenv('SLURM_ARRAY_TASK_ID')));
if ~isfinite(taskId)
    taskId = 1;
end
taskId = round(taskId);
setupFile = getenv('FPP_SETUP_FILE');
outputRoot = getenv('FPP_OUTPUT_ROOT');
if isempty(setupFile) || ~isfile(setupFile)
    error('FPP:SetupFile', 'FPP_SETUP_FILE must identify global_bifurcation_setup.mat.');
end
if isempty(outputRoot)
    error('FPP:OutputRoot', 'FPP_OUTPUT_ROOT is required.');
end
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

loaded = load(setupFile, 'setup');
setup = loaded.setup;
caseTable = local_case_table(setup);
if taskId < 1 || taskId > height(caseTable)
    error('FPP:TaskId', 'Task id %d is outside 1:%d.', taskId, height(caseTable));
end
spec = caseTable(taskId, :);

context = setup.Context;
fixed = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
if numel(fixed) ~= 3*n
    error('FPP:MapSize', 'Fixed point and map size are inconsistent.');
end

if spec.pathway == "L6"
    gamma6 = 1-spec.weight;
    gammaI = 1;
    jacobian = setup.Pathway.JRest + gamma6*setup.Pathway.J6 + setup.Pathway.JI;
    phi = @(state) l6ns_phi(state, spec.weight, context, [1 1], 1, 'fpp');
    criticalDirection = real(setup.L6.Right(:));
    criticalWeight = setup.L6.W;
else
    gamma6 = 1;
    gammaI = 1-spec.weight;
    jacobian = setup.Pathway.JRest + setup.Pathway.J6 + gammaI*setup.Pathway.JI;
    phi = @(state) l6ns_phi(state, 0, context, [1 1], gammaI, 'fpp');
    criticalDirection = real(setup.Inhibition.Right(:));
    criticalWeight = setup.Inhibition.W;
end

[maxRealLambda, leadingLambda] = local_max_real(jacobian);
fixedResidual = norm(phi(fixed)-fixed)/max(norm(fixed), eps);
if fixedResidual > 1e-10
    error('FPP:FixedPoint', ...
        'Fixed-point-preserving residual %.3e is larger than tolerance.', fixedResidual);
end

criticalDirection = criticalDirection/norm(criticalDirection);
[~, orientationIndex] = max(abs(criticalDirection));
if criticalDirection(orientationIndex) < 0
    criticalDirection = -criticalDirection;
end
requestedInitialHC = 0.01;
unitHC = local_hcnorm(fixed+criticalDirection, fixed, n, context.CWeight);
amplitude = requestedInitialHC/max(unitHC, eps);
negative = criticalDirection < 0;
if any(negative)
    nonnegativeLimit = 0.5*min(fixed(negative)./(-criticalDirection(negative)));
    amplitude = min(amplitude, nonnegativeLimit);
end
if ~isfinite(amplitude) || amplitude <= 0
    error('FPP:Perturbation', 'Could not construct a positive initial perturbation.');
end
initialState = fixed + amplitude*criticalDirection;
initialHC = local_hcnorm(initialState, fixed, n, context.CWeight);

tauMs = setup.Config.TauMs;
if spec.regime == "stable-control"
    finalTimeMs = 1500;
else
    growth = maxRealLambda-1;
    if growth <= 0
        error('FPP:StabilitySide', ...
            'Expected an unstable case, but max Re lambda is %.9f.', maxRealLambda);
    end
    finalTimeMs = min(8000, max(300, 18*tauMs/growth));
end
snapshotTimes = linspace(finalTimeMs/10, finalTimeMs, 10);

blowupThresholdHz = 500;
rhs = @(~,state) (phi(state)-state)/tauMs;
maxStepMs = max(0.25, min(10, finalTimeMs/100));
options = odeset('RelTol', 2e-5, 'AbsTol', 1e-7, ...
    'MaxStep', maxStepMs, 'Events', @local_blowup_event);
fprintf(['FPP task %d/%d: %s %s, w=%+.9f, gamma6=%.9f, gammaI=%.9f, ' ...
    'maxRe=%.9f, T=%.3f ms.\n'], taskId, height(caseTable), spec.pathway, ...
    spec.regime, spec.weight, gamma6, gammaI, maxRealLambda, finalTimeMs);
solution = ode45(rhs, [0 finalTimeMs], initialState, options);

sampleTimes = [0 snapshotTimes];
sampleStates = nan(numel(fixed), numel(sampleTimes));
valid = sampleTimes <= solution.x(end) + 10*eps(max(1,solution.x(end)));
if any(valid)
    sampleStates(:,valid) = deval(solution, sampleTimes(valid));
end
thresholdCrossed = solution.x(end) < finalTimeMs-1e-8*max(1,finalTimeMs);

hc = nan(size(sampleTimes));
relativeResidual = nan(size(sampleTimes));
minimumRate = nan(size(sampleTimes));
maximumRate = nan(size(sampleTimes));
for index = 1:numel(sampleTimes)
    state = sampleStates(:,index);
    if all(isfinite(state))
        hc(index) = local_hcnorm(state, fixed, n, context.CWeight);
        relativeResidual(index) = norm(phi(state)-state)/max(norm(state),eps);
        minimumRate(index) = min(state);
        maximumRate(index) = max(state);
    end
end

caseTag = sprintf('%02d_%s_%s_w%s', taskId, lower(char(spec.pathway)), ...
    strrep(char(spec.regime), '-', '_'), local_number_tag(spec.weight));
caseDir = fullfile(outputRoot, caseTag);
if ~exist(caseDir, 'dir')
    mkdir(caseDir);
end

snapshotTable = table(repmat(taskId,numel(sampleTimes),1), ...
    repmat(spec.pathway,numel(sampleTimes),1), ...
    repmat(spec.regime,numel(sampleTimes),1), ...
    repmat(spec.weight,numel(sampleTimes),1), ...
    (0:numel(sampleTimes)-1)', sampleTimes(:), hc(:), relativeResidual(:), ...
    minimumRate(:), maximumRate(:), valid(:), ...
    'VariableNames', {'taskId','pathway','regime','weight','snapshotIndex', ...
    'timeMs','HCnormFromBaseline','relativeMapResidual','minimumRateHz', ...
    'maximumRateHz','available'});
writetable(snapshotTable, fullfile(caseDir, 'snapshot_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

finalAvailable = find(valid, 1, 'last');
caseSummary = table(taskId, spec.pathway, spec.regime, spec.weight, ...
    criticalWeight, gamma6, gammaI, real(leadingLambda), imag(leadingLambda), ...
    maxRealLambda, tauMs, finalTimeMs, solution.x(end), requestedInitialHC, ...
    initialHC, fixedResidual, thresholdCrossed, blowupThresholdHz, ...
    hc(finalAvailable), minimumRate(finalAvailable), maximumRate(finalAvailable), ...
    'VariableNames', {'taskId','pathway','regime','weight','criticalWeight', ...
    'gamma6','gammaI','leadingLambdaReal','leadingLambdaImag','maxRealLambda', ...
    'tauMs','requestedFinalTimeMs','actualFinalTimeMs','requestedInitialHCnorm', ...
    'initialHCnorm','fixedPointResidual','blowupThresholdCrossed', ...
    'blowupThresholdHz','lastHCnorm','lastMinimumRateHz','lastMaximumRateHz'});
writetable(caseSummary, fullfile(caseDir, 'case_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

metadata = struct('TaskId',taskId,'Spec',spec,'CriticalWeight',criticalWeight, ...
    'Gamma6',gamma6,'GammaI',gammaI,'LeadingLambda',leadingLambda, ...
    'MaxRealLambda',maxRealLambda,'TauMs',tauMs,'FinalTimeMs',finalTimeMs, ...
    'InitialAmplitude',amplitude,'InitialHCnorm',initialHC, ...
    'FixedPointResidual',fixedResidual,'ThresholdCrossed',thresholdCrossed, ...
    'BlowupThresholdHz',blowupThresholdHz);
save(fullfile(caseDir, 'trajectory_result.mat'), 'metadata', 'sampleTimes', ...
    'sampleStates', 'hc', 'relativeResidual', 'minimumRate', 'maximumRate', ...
    '-v7.3');

if ~strcmp(getenv('FPP_SKIP_PLOT'),'1')
    if spec.regime == "stable-control"
        local_plot_stable(fixed, sampleStates(:,finalAvailable), mapSize, context.CWeight, ...
            caseSummary, hc(finalAvailable), caseDir);
    else
        local_plot_trajectory(fixed, sampleStates(:,2:end), snapshotTimes, hc(2:end), ...
            valid(2:end), mapSize, context.CWeight, caseSummary, caseDir);
    end
end
fprintf('Completed %s: last HC=%.6g, rate range [%.6g, %.6g] Hz, threshold=%d.\n', ...
    caseTag, hc(finalAvailable), minimumRate(finalAvailable), ...
    maximumRate(finalAvailable), thresholdCrossed);

    function [value,isterminal,direction] = local_blowup_event(~,state)
        value = blowupThresholdHz-max(abs(state));
        isterminal = 1;
        direction = -1;
    end
end

function cases = local_case_table(setup)
cases = table((1:10)', ...
    [repmat("L6",5,1); repmat("Inhibition",5,1)], ...
    ["stable-control"; repmat("unstable",4,1); ...
     "stable-control"; repmat("unstable",4,1)], ...
    [-0.15; -0.20; -0.225; -0.25; -0.30; ...
      0.125; 0.175; 0.20; 0.225; 0.25], ...
    'VariableNames', {'taskId','pathway','regime','weight'});
if ~(setup.L6.W < -0.15 && setup.L6.W > -0.20)
    error('FPP:L6Boundary', 'Unexpected L6 critical weight %.12g.', setup.L6.W);
end
if ~(setup.Inhibition.W > 0.125 && setup.Inhibition.W < 0.175)
    error('FPP:IBoundary', ...
        'Unexpected inhibition critical weight %.12g.', setup.Inhibition.W);
end
end

function [maximumReal, leadingLambda] = local_max_real(jacobian)
options = struct('tol',1e-9,'maxit',1600,'p',80,'isreal',true);
try
    values = eigs(jacobian, 8, 'largestreal', options);
catch
    values = eigs(jacobian, 8, 'lr', options);
end
[maximumReal,index] = max(real(values));
leadingLambda = values(index);
end

function value = local_hcnorm(state, baseline, n, wC)
state = state(:);
baseline = baseline(:);
dS = state(1:n)-baseline(1:n);
dC = state(n+(1:n))-baseline(n+(1:n));
dI = state(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-wC)*dS+wC*dC;
value = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end

function maps = local_state_maps(state, mapSize, wC)
n = prod(mapSize);
maps = cell(4,1);
maps{1} = reshape(state(1:n),mapSize);
maps{2} = reshape(state(n+(1:n)),mapSize);
maps{3} = reshape(state(2*n+(1:n)),mapSize);
maps{4} = (1-wC)*maps{1}+wC*maps{2};
end

function local_plot_baseline_or_state(axisHandle, map, limits)
imagesc(axisHandle, map);
axis(axisHandle, 'image');
axis(axisHandle, 'off');
clim(axisHandle, limits);
end

function local_plot_stable(fixed, finalState, mapSize, wC, summary, finalHC, outputDir)
populationNames = {'S','C','I','E'};
states = {fixed, finalState};
allMaps = cell(4,2);
for column = 1:2
    allMaps(:,column) = local_state_maps(states{column},mapSize,wC);
end
limits = local_population_limits(allMaps);
fig = figure('Visible','off','Color','w','Position',[40 40 1150 1700]);
layout = tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');
for row = 1:4
    for column = 1:2
        ax = nexttile(layout,(row-1)*2+column);
        local_plot_baseline_or_state(ax,allMaps{row,column},limits(row,:));
        if row == 1
            if column == 1
                title(ax,'Baseline fixed point, HC=0','FontSize',12);
            else
                title(ax,sprintf('Final t=%.1f ms, HC=%.5g', ...
                    summary.actualFinalTimeMs,finalHC),'FontSize',12);
            end
        end
        if column == 1
            text(ax,-0.08,0.5,populationNames{row},'Units','normalized', ...
                'HorizontalAlignment','right','FontWeight','bold','FontSize',13);
        else
            cb = colorbar(ax,'eastoutside');
            cb.FontSize = 11;
        end
    end
end
colormap(fig,parula(256));
title(layout,sprintf(['%s stable-side FPP ODE control: w=%+.3f, w_c=%+.6f, ' ...
    'max Re lambda=%.6f, tau=%.4f ms'], summary.pathway,summary.weight, ...
    summary.criticalWeight,summary.maxRealLambda,summary.tauMs), ...
    'Interpreter','none','FontWeight','bold','FontSize',14);
exportgraphics(fig,fullfile(outputDir,'stable_side_baseline_vs_final.pdf'), ...
    'ContentType','image','Resolution',220);
close(fig);
end

function local_plot_trajectory(fixed, snapshotStates, snapshotTimes, hc, available, ...
        mapSize, wC, summary, outputDir)
populationNames = {'S','C','I','E'};
stateCount = 1+numel(snapshotTimes);
allMaps = cell(4,stateCount);
allMaps(:,1) = local_state_maps(fixed,mapSize,wC);
for column = 2:stateCount
    if available(column-1) && all(isfinite(snapshotStates(:,column-1)))
        allMaps(:,column) = local_state_maps(snapshotStates(:,column-1),mapSize,wC);
    else
        allMaps(:,column) = repmat({nan(mapSize)},4,1);
    end
end
limits = local_population_limits(allMaps);
fig = figure('Visible','off','Color','w','Position',[20 20 4600 1750]);
layout = tiledlayout(fig,4,stateCount,'TileSpacing','compact','Padding','compact');
for row = 1:4
    for column = 1:stateCount
        ax = nexttile(layout,(row-1)*stateCount+column);
        local_plot_baseline_or_state(ax,allMaps{row,column},limits(row,:));
        if row == 1
            if column == 1
                title(ax,'Baseline f*, HC=0','FontSize',16);
            elseif available(column-1)
                title(ax,sprintf('t=%.1f ms\nHC=%.4g', ...
                    snapshotTimes(column-1),hc(column-1)),'FontSize',16);
            else
                title(ax,sprintf('t=%.1f ms\nterminated', ...
                    snapshotTimes(column-1)),'FontSize',16);
            end
        end
        if column == 1
            text(ax,-0.08,0.5,populationNames{row},'Units','normalized', ...
                'HorizontalAlignment','right','FontWeight','bold','FontSize',18);
        elseif column == stateCount
            cb = colorbar(ax,'eastoutside');
            cb.FontSize = 14;
        end
    end
end
colormap(fig,parula(256));
title(layout,sprintf(['%s unstable FPP ODE: w=%+.3f, w_c=%+.6f, ' ...
    'max Re lambda=%.6f, tau=%.4f ms, threshold crossed=%d'], ...
    summary.pathway,summary.weight,summary.criticalWeight, ...
    summary.maxRealLambda,summary.tauMs,summary.blowupThresholdCrossed), ...
    'Interpreter','none','FontWeight','bold','FontSize',22);
exportgraphics(fig,fullfile(outputDir,'baseline_and_10_ode_snapshots.pdf'), ...
    'ContentType','image','Resolution',180);
close(fig);
end

function limits = local_population_limits(allMaps)
limits = nan(4,2);
for row = 1:4
    values = [];
    for column = 1:size(allMaps,2)
        map = allMaps{row,column};
        values = [values; map(isfinite(map))]; %#ok<AGROW>
    end
    if isempty(values)
        limits(row,:) = [0 1];
    else
        lower = min(values);
        upper = max(values);
        if upper <= lower
            padding = max(1e-6,abs(lower)*1e-6);
            lower = lower-padding;
            upper = upper+padding;
        end
        limits(row,:) = [lower upper];
    end
end
end

function value = local_env_number(name, defaultValue)
raw = getenv(name);
if isempty(raw)
    value = defaultValue;
else
    value = str2double(raw);
end
end

function tag = local_number_tag(value)
tag = sprintf('%+.3f',value);
tag = strrep(tag,'+','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'.','p');
end
