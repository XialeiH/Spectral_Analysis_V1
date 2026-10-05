function run_real_moved_iteration_task()
% Run one real whole-pathway regular-iteration trajectory.

taskId = local_env_number('REAL_MOVED_TASK_ID', ...
    str2double(getenv('SLURM_ARRAY_TASK_ID')));
if ~isfinite(taskId)
    taskId = 1;
end
taskId = round(taskId);
setupFile = getenv('REAL_MOVED_SETUP_FILE');
outputRoot = getenv('REAL_MOVED_OUTPUT_ROOT');
if isempty(setupFile) || ~isfile(setupFile)
    error('RealMoved:SetupFile', ...
        'REAL_MOVED_SETUP_FILE must identify global_bifurcation_setup.mat.');
end
if isempty(outputRoot)
    error('RealMoved:OutputRoot', 'REAL_MOVED_OUTPUT_ROOT is required.');
end
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

loaded = load(setupFile, 'setup');
setup = loaded.setup;
caseTable = local_case_table();
if taskId < 1 || taskId > height(caseTable)
    error('RealMoved:TaskId', 'Task id %d is outside 1:%d.', ...
        taskId, height(caseTable));
end
spec = caseTable(taskId, :);

context = setup.Context;
baseline = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
if numel(baseline) ~= 3*n
    error('RealMoved:MapSize', 'Fixed point and map size are inconsistent.');
end

gain6 = 1;
gainI = 1;
if spec.pathway == "L6"
    gain6 = spec.pathwayGain;
else
    gainI = spec.pathwayGain;
end
phi = @(state) l6ns_phi(state, 1-gain6, context, [1 1], gainI, 'true');

fprintf('Real moved-FP iteration task %d/%d: %s beta=%+.4f, gain=%.4f.\n', ...
    taskId, height(caseTable), spec.pathway, spec.beta, spec.pathwayGain);
movedResult = real_tuning_fixed_point(phi, baseline, context.RelaxationP);
if ~movedResult.Converged
    error('RealMoved:FixedPoint', ...
        'Moved-point iteration did not converge: %s, residual %.3e.', ...
        movedResult.Termination, movedResult.Residual);
end
moved = movedResult.State(:);
movedResidual = norm(phi(moved)-moved)/max(norm(moved), eps);

jacobian = real_tuning_true_jacobian(moved, context, gain6, gainI);
[maxRealLambda, leadingLambda] = local_max_real(jacobian);
if maxRealLambda >= 1
    error('RealMoved:UnstableMovedPoint', ...
        'The moved fixed point has max Re lambda %.9f >= 1.', maxRealLambda);
end

[sampleIterations,sampleStates] = local_iteration_snapshots(phi, baseline, ...
    context.RelaxationP, movedResult.Iterations);
finalRelativeError = norm(sampleStates(:,end)-moved)/max(norm(moved),eps);
if finalRelativeError > 1e-12
    error('RealMoved:IterationReplay', ...
        'Replayed iteration differs from the moved fixed point by %.3e.', ...
        finalRelativeError);
end

hcFromBaseline = nan(size(sampleIterations));
hcToMoved = nan(size(sampleIterations));
relativeMapResidual = nan(size(sampleIterations));
relativeDistanceToMoved = nan(size(sampleIterations));
minimumRate = nan(size(sampleIterations));
maximumRate = nan(size(sampleIterations));
for index = 1:numel(sampleIterations)
    state = sampleStates(:,index);
    hcFromBaseline(index) = local_hcnorm(state, baseline, n, context.CWeight);
    hcToMoved(index) = local_hcnorm(state, moved, n, context.CWeight);
    relativeMapResidual(index) = norm(phi(state)-state)/max(norm(state),eps);
    relativeDistanceToMoved(index) = norm(state-moved)/max(norm(moved),eps);
    minimumRate(index) = min(state);
    maximumRate(index) = max(state);
end

caseTag = sprintf('%02d_%s_real_beta_%s', taskId, ...
    lower(char(spec.pathway)), local_number_tag(spec.beta));
caseDir = fullfile(outputRoot, caseTag);
if ~exist(caseDir, 'dir')
    mkdir(caseDir);
end

snapshotTable = table(repmat(taskId,numel(sampleIterations),1), ...
    repmat(spec.pathway,numel(sampleIterations),1), ...
    repmat(spec.beta,numel(sampleIterations),1), ...
    repmat(spec.pathwayGain,numel(sampleIterations),1), ...
    (0:numel(sampleIterations)-1)', sampleIterations(:), ...
    hcFromBaseline(:), hcToMoved(:), relativeMapResidual(:), ...
    relativeDistanceToMoved(:), minimumRate(:), maximumRate(:), ...
    'VariableNames', {'taskId','pathway','beta','pathwayGain', ...
    'snapshotIndex','iteration','HCnormFromOriginalBaseline', ...
    'HCnormToMovedFixedPoint','relativeMapResidual', ...
    'relativeDistanceToMovedFixedPoint','minimumRateHz','maximumRateHz'});
writetable(snapshotTable, fullfile(caseDir, 'snapshot_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

movedHC = local_hcnorm(moved, baseline, n, context.CWeight);
caseSummary = table(taskId, spec.pathway, spec.beta, spec.pathwayGain, ...
    -spec.beta, gain6, gainI, context.RelaxationP, ...
    movedResult.Iterations, movedResidual, movedHC, ...
    movedResult.MinimumRateHz, movedResult.MaximumRateHz, ...
    real(leadingLambda), imag(leadingLambda), maxRealLambda, ...
    finalRelativeError, ...
    'VariableNames', {'taskId','pathway','beta','pathwayGain','matchedFppW', ...
    'gain6','gainI','relaxationP','fixedPointIterations', ...
    'fixedPointResidual','movedHCnormFromOriginalBaseline', ...
    'movedMinimumRateHz','movedMaximumRateHz','leadingLambdaReal', ...
    'leadingLambdaImag','maxRealLambda','iterationReplayRelativeError'});
writetable(caseSummary, fullfile(caseDir, 'case_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

metadata = struct('TaskId',taskId,'Spec',spec,'Gain6',gain6, ...
    'GainI',gainI,'MatchedFppW',-spec.beta,'LeadingLambda',leadingLambda, ...
    'MaxRealLambda',maxRealLambda,'RelaxationP',context.RelaxationP, ...
    'MovedFixedPointResult',movedResult, ...
    'MovedFixedPointResidual',movedResidual, ...
    'IterationReplayRelativeError',finalRelativeError);
save(fullfile(caseDir, 'iteration_result.mat'), 'metadata', 'moved', ...
    'sampleIterations', 'sampleStates', 'hcFromBaseline', 'hcToMoved', ...
    'relativeMapResidual', 'relativeDistanceToMoved', ...
    'minimumRate', 'maximumRate', '-v7.3');

if spec.summaryLayout == "pair"
    local_plot_pair(baseline, moved, mapSize, context.CWeight, ...
        caseSummary, caseDir);
else
    local_plot_trajectory(baseline, sampleStates(:,2:end), ...
        sampleIterations(2:end), hcFromBaseline(2:end), hcToMoved(2:end), ...
        mapSize, context.CWeight, caseSummary, caseDir);
end
fprintf(['Completed %s: iterations=%d, moved HC=%.6g, residual=%.3e, ' ...
    'replay error=%.3e.\n'], caseTag, movedResult.Iterations, movedHC, ...
    movedResidual, finalRelativeError);
end

function [sampleIterations,sampleStates] = local_iteration_snapshots( ...
        phi, baseline, relaxationP, totalIterations)
snapshotIterations = round(linspace(1,totalIterations,10));
if numel(unique(snapshotIterations)) ~= 10
    error('RealMoved:Snapshots', ...
        'Need at least ten distinct regular-iteration snapshots.');
end
sampleIterations = [0 snapshotIterations];
sampleStates = nan(numel(baseline),numel(sampleIterations));
sampleStates(:,1) = baseline;
state = baseline;
response = phi(state);
snapshotColumn = 2;
for iteration = 1:totalIterations
    state = (1-relaxationP)*state + relaxationP*response;
    response = phi(state);
    if iteration == snapshotIterations(snapshotColumn-1)
        sampleStates(:,snapshotColumn) = state;
        snapshotColumn = snapshotColumn+1;
        if snapshotColumn > numel(sampleIterations)
            break
        end
    end
end
if any(~isfinite(sampleStates),'all')
    error('RealMoved:Snapshots', 'Iteration replay produced missing snapshots.');
end
end

function cases = local_case_table()
cases = table((1:16)', ...
    [repmat("L6",5,1); repmat("Inhibition",11,1)], ...
    [0.15; 0.20; 0.225; 0.25; 0.30; ...
    -0.125; -0.175; -0.20; -0.225; -0.25; ...
    -0.05; -0.10; -0.15; -0.075; -0.085; -0.08], ...
    [1.15; 1.20; 1.225; 1.25; 1.30; ...
    0.875; 0.825; 0.80; 0.775; 0.75; ...
    0.95; 0.90; 0.85; 0.925; 0.915; 0.92], ...
    ["pair"; repmat("trajectory",4,1); ...
    "pair"; repmat("trajectory",4,1); repmat("pair",6,1)], ...
    'VariableNames', {'taskId','pathway','beta','pathwayGain','summaryLayout'});
end

function [maximumReal, leadingLambda] = local_max_real(jacobian)
options = struct('tol',1e-9,'maxit',1800,'p',80,'isreal',true,'disp',0);
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

function local_plot_state(axisHandle, map, limits)
imagesc(axisHandle, map);
axis(axisHandle, 'image');
axis(axisHandle, 'off');
clim(axisHandle, limits);
end

function local_plot_pair(baseline, moved, mapSize, wC, summary, outputDir)
populationNames = {'S','C','I','E'};
states = {baseline, moved};
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
        local_plot_state(ax,allMaps{row,column},limits(row,:));
        if row == 1
            if column == 1
                title(ax,'Original baseline fixed point, HC=0','FontSize',12);
            else
                title(ax,sprintf(['Moved fixed point after %d iterations\n' ...
                    'HC(base)=%.5g, residual=%.2e'], ...
                    summary.fixedPointIterations, ...
                    summary.movedHCnormFromOriginalBaseline, ...
                    summary.fixedPointResidual),'FontSize',12);
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
title(layout,sprintf(['%s real moved-fixed-point iteration: beta=%+.3f, ' ...
    'gain=%.3f, max Re lambda=%.6f, p=%.2f'], ...
    summary.pathway,summary.beta,summary.pathwayGain,summary.maxRealLambda, ...
    summary.relaxationP),'Interpreter','none','FontWeight','bold','FontSize',14);
exportgraphics(fig,fullfile(outputDir,'baseline_vs_moved_fixed_point.pdf'), ...
    'ContentType','image','Resolution',220);
close(fig);
end

function local_plot_trajectory(baseline, snapshotStates, snapshotIterations, ...
        hcBase, hcMoved, mapSize, wC, summary, outputDir)
populationNames = {'S','C','I','E'};
stateCount = 1+numel(snapshotIterations);
allMaps = cell(4,stateCount);
allMaps(:,1) = local_state_maps(baseline,mapSize,wC);
for column = 2:stateCount
    allMaps(:,column) = local_state_maps(snapshotStates(:,column-1),mapSize,wC);
end
limits = local_population_limits(allMaps);
fig = figure('Visible','off','Color','w','Position',[20 20 4600 1750]);
layout = tiledlayout(fig,4,stateCount,'TileSpacing','compact','Padding','compact');
for row = 1:4
    for column = 1:stateCount
        ax = nexttile(layout,(row-1)*stateCount+column);
        local_plot_state(ax,allMaps{row,column},limits(row,:));
        if row == 1
            if column == 1
                title(ax,'Original baseline, HC=0','FontSize',16);
            else
                title(ax,sprintf(['iteration %d\nHC(base)=%.4g\n' ...
                    'HC(target)=%.3g'], snapshotIterations(column-1), ...
                    hcBase(column-1),hcMoved(column-1)),'FontSize',14);
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
title(layout,sprintf(['%s real moved-fixed-point iteration: beta=%+.3f, ' ...
    'gain=%.3f, max Re lambda=%.6f, p=%.2f, final residual=%.2e'], ...
    summary.pathway,summary.beta,summary.pathwayGain,summary.maxRealLambda, ...
    summary.relaxationP,summary.fixedPointResidual), ...
    'Interpreter','none','FontWeight','bold','FontSize',22);
exportgraphics(fig,fullfile(outputDir, ...
    'baseline_and_10_iteration_snapshots.pdf'), ...
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
