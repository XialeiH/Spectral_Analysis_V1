function summary = run_inhibition_multiangle_top4_atlas_task(taskIndex, sourceRoot, sweepRoot, outputRoot)
% Plot top-four positive-edge modes for inhibition components and selected J(w_I).

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('INHATLAS_SOURCE_ROOT');
end
if nargin < 3 || isempty(sweepRoot)
    sweepRoot = getenv('INHATLAS_SWEEP_ROOT');
end
if nargin < 4 || isempty(outputRoot)
    outputRoot = getenv('INHATLAS_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrast = 100;
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > numel(angles)
    error('InhibitionAtlas:TaskIndex', 'Task index must be an integer from 1 to 4.');
end
if isempty(sourceRoot) || ~isfolder(sourceRoot)
    error('InhibitionAtlas:SourceRoot', ...
        'INHATLAS_SOURCE_ROOT must identify the audited cell-type endpoint root.');
end
if isempty(sweepRoot) || ~isfolder(sweepRoot)
    error('InhibitionAtlas:SweepRoot', ...
        'INHATLAS_SWEEP_ROOT must identify the archived inhibition sweep results.');
end
if isempty(outputRoot)
    error('InhibitionAtlas:OutputRoot', 'INHATLAS_OUTPUT_ROOT is required.');
end

angleValue = angles(taskIndex);
[jacobian0, jacobian1] = local_endpoint_jacobians(sourceRoot, angleValue, contrast);
jBase = sparse(jacobian1);
jI = sparse(jacobian0 - jacobian1);
clear jacobian0 jacobian1

boundary = local_boundary(sweepRoot, angleValue, contrast);
[stableW, stableMaxReal, stableMargin] = local_highly_stable_weight( ...
    sweepRoot, angleValue, contrast);
if boundary.lowMaxReal >= 1 || boundary.highMaxReal < 1
    error('InhibitionAtlas:BoundarySides', ...
        'The stored bracket does not place lowW below and highW above lambda=1.');
end
if stableMargin <= 0
    error('InhibitionAtlas:StableWeight', ...
        'The selected highly stable inhibition weight is not stable.');
end

cfg = l6ns_config();
cfg.EigsTolerance = 1e-10;
cfg.EigsMaxIterations = 2200;
cfg.EigsSubspaceDimension = 80;
cfg.RandomSeed = 2200 + taskIndex;
rng(cfg.RandomSeed, 'twister');

figureRoot = fullfile(outputRoot, 'figures');
dataRoot = fullfile(outputRoot, 'data');
if ~exist(figureRoot, 'dir'); mkdir(figureRoot); end
if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end

summaryRows = {};
modeData = struct();

componentSpecs = { ...
    'J_base_I', jBase, 1, ...
        'A_I = J(w_I=1) = J_rest+J_6, inhibition-independent base'; ...
    'J_I', jI, NaN, ...
        'B_I = J(w_I=0)-J(w_I=1), w_I-independent inhibition component'};
for componentIndex = 1:size(componentSpecs, 1)
    focus = componentSpecs{componentIndex, 1};
    matrix = componentSpecs{componentIndex, 2};
    displayedW = componentSpecs{componentIndex, 3};
    definition = componentSpecs{componentIndex, 4};
    modes = l6ns_eigenpairs(matrix, 4, 'largestreal', cfg);
    if strcmp(focus, 'J_I') && max(abs(modes.Lambda)) < 1e-6
        definition = [definition, ...
            ' | near-zero nullspace edge (degenerate; basis-dependent)'];
    end
    stateLabel = 'Component_Positive_Edge';
    titleLines = { ...
        sprintf('%s positive-edge top 4 | angle %.2f deg | contrast %d', ...
        focus, angleValue, contrast), ...
        sprintf('%s | displayed w_I=%s | max Re(lambda)=%.9f', ...
        definition, local_display_number(displayedW), real(modes.Lambda(1)))};
    stem = sprintf('Top4_%s_Angle_%sdeg_Contrast_%d', ...
        focus, local_angle_tag(angleValue), contrast);
    local_plot_mode_atlas(modes, angleValue, contrast, focus, displayedW, ...
        titleLines, figureRoot, stem);
    [rows, stored] = local_mode_records(modes, angleValue, contrast, focus, ...
        stateLabel, displayedW, NaN, stem);
    summaryRows = [summaryRows; rows]; %#ok<AGROW>
    modeData.(focus) = stored;
end

fullW = [0, 1, 0.5, stableW, boundary.lowW, boundary.highW];
stateLabels = {'Baseline_Dynamic_Inhibition_wI0', 'Frozen_Inhibition_wI1', ...
    'Partial_Freeze_wI0p5', 'Highly_Stable_Stronger_Inhibition', ...
    'Boundary_Just_Stable', 'Boundary_Just_Unstable'};
for stateIndex = 1:numel(fullW)
    wValue = fullW(stateIndex);
    gammaI = 1 - wValue;
    matrix = jBase + gammaI * jI;
    modes = l6ns_eigenpairs(matrix, 4, 'largestreal', cfg);
    maxReal = real(modes.Lambda(1));
    margin = 1 - maxReal;
    status = 'stable';
    if maxReal >= 1; status = 'unstable'; end
    focus = 'Full_JwI';
    stateLabel = stateLabels{stateIndex};
    titleLines = { ...
        sprintf('Full J(w_I) positive-edge top 4 | angle %.2f deg | contrast %d | %s', ...
        angleValue, contrast, strrep(stateLabel, '_', ' ')), ...
        sprintf(['J(w_I)=J_base_I+(1-w_I)J_I | w_I=%+.9f | gamma_I=%.9f | ' ...
        'max Re(lambda)=%.9f | margin=%.9f (%s)'], ...
        wValue, gammaI, maxReal, margin, status)};
    stem = sprintf('Top4_Full_JwI_Angle_%sdeg_Contrast_%d_%s_wI%s', ...
        local_angle_tag(angleValue), contrast, stateLabel, local_number_tag(wValue, 9));
    local_plot_mode_atlas(modes, angleValue, contrast, focus, wValue, ...
        titleLines, figureRoot, stem);
    [rows, stored] = local_mode_records(modes, angleValue, contrast, focus, ...
        stateLabel, wValue, gammaI, stem);
    summaryRows = [summaryRows; rows]; %#ok<AGROW>
    modeData.(stateLabel) = stored;
end

summary = cell2table(summaryRows, 'VariableNames', { ...
    'angle','contrast','jacobianFocus','state','wI','gammaI','mode', ...
    'lambdaReal','lambdaImag','conditionNumber','rightResidual','leftResidual', ...
    'S_energyFraction','C_energyFraction','I_energyFraction','figureStem'});
summaryFile = fullfile(dataRoot, sprintf( ...
    'top4_inhibition_mode_summary_angle%s_contrast%d.tsv', ...
    local_angle_tag(angleValue), contrast));
writetable(summary, summaryFile, 'FileType', 'text', 'Delimiter', '\t');

metadata = struct( ...
    'JacobianDefinition', 'J(w_I)=J_base_I+(1-w_I)*J_I', ...
    'JBaseDefinition', 'J_base_I=J(w_I=1)=J_rest+J_6', ...
    'JIDefinition', 'J_I=J(w_I=0)-J(w_I=1)', ...
    'JIPositiveEdgeInterpretation', ...
        ['The positive edge of J_I is a degenerate near-zero nullspace; ' ...
         'the four displayed vectors are one numerical basis, not unique modes.'], ...
    'Angle', angleValue, 'Contrast', contrast, ...
    'BoundaryWI', boundary.wBoundary, ...
    'BoundaryLowW', boundary.lowW, ...
    'BoundaryHighW', boundary.highW, ...
    'HighlyStableWI', stableW, ...
    'HighlyStableMaxReal', stableMaxReal, ...
    'HighlyStableMargin', stableMargin, ...
    'SavedFullMatrices', false);
save(fullfile(dataRoot, sprintf( ...
    'top4_inhibition_mode_data_angle%s_contrast%d.mat', ...
    local_angle_tag(angleValue), contrast)), ...
    'summary', 'metadata', 'modeData', '-v7.3');
fprintf('Saved eight inhibition top-four atlases for angle %.2f to %s.\n', ...
    angleValue, outputRoot);
end

function [jacobian0, jacobian1] = local_endpoint_jacobians(sourceRoot, angleValue, contrast)
file0 = local_endpoint_file(sourceRoot, angleValue, contrast, 0);
file1 = local_endpoint_file(sourceRoot, angleValue, contrast, 1);
loaded0 = load(file0, 'Section4');
loaded1 = load(file1, 'Section4');
jacobian0 = loaded0.Section4.A;
jacobian1 = loaded1.Section4.A;
if ~isequal(size(jacobian0), size(jacobian1)) || size(jacobian0,1) ~= size(jacobian0,2)
    error('InhibitionAtlas:JacobianSize', ...
        'Endpoint Jacobians have incompatible dimensions.');
end
end

function file = local_endpoint_file(sourceRoot, angleValue, contrast, wValue)
tag = sprintf('L6eqWS0p00_C0p00_I%.2f', wValue);
tag = strrep(strrep(tag, '.', 'p'), '-', 'm');
file = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_%s_contr%d_angle_%.2f.mat', ...
    tag, contrast, angleValue));
if ~isfile(file)
    error('InhibitionAtlas:MissingEndpoint', 'Missing endpoint file %s.', file);
end
end

function boundary = local_boundary(sweepRoot, angleValue, contrast)
file = fullfile(sweepRoot, 'boundary_data', sprintf( ...
    'inhibition_boundary_angle%.2f_contrast%d.tsv', angleValue, contrast));
if ~isfile(file)
    error('InhibitionAtlas:MissingBoundary', 'Missing exact boundary file %s.', file);
end
tableIn = readtable(file, 'FileType', 'text', 'Delimiter', '\t');
if height(tableIn) ~= 1
    error('InhibitionAtlas:BoundaryRows', ...
        'Expected exactly one boundary row in %s.', file);
end
boundary = table2struct(tableIn);
end

function [wValue, maxReal, margin] = local_highly_stable_weight(sweepRoot, angleValue, contrast)
file = fullfile(sweepRoot, 'sweep_data', sprintf( ...
    'inhibition_weight_sweep_angle%.2f_contrast%d.tsv', angleValue, contrast));
if ~isfile(file)
    error('InhibitionAtlas:MissingSweep', 'Missing archived sweep file %s.', file);
end
tableIn = readtable(file, 'FileType', 'text', 'Delimiter', '\t');
stable = tableIn.stabilityMargin > 0;
if ~any(stable)
    error('InhibitionAtlas:NoStableSweepPoint', ...
        'No stable point was found in %s.', file);
end
stableRows = find(stable);
[margin, localIndex] = max(tableIn.stabilityMargin(stable));
row = stableRows(localIndex);
wValue = tableIn.wI(row);
maxReal = tableIn.maxRealLambda(row);
end

function local_plot_mode_atlas(modes, angleValue, contrast, focus, wValue, ...
        titleLines, outputRoot, stem)
n = size(modes.Right, 1) / 3;
side = round(sqrt(n));
if n ~= round(n) || side * side ~= n
    error('InhibitionAtlas:MapSize', ...
        'The three-population state does not form square maps.');
end
mapSize = [side side];
populationLabels = {'S','C','I','E = 0.6923 S + 0.3077 C'};

figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [40 40 2200 1850]);
layout = tiledlayout(figureHandle, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
colormap(figureHandle, parula(256));
for modeIndex = 1:4
    vector = modes.Right(:, modeIndex);
    s = vector(1:n);
    c = vector(n + (1:n));
    i = vector(2*n + (1:n));
    e = 0.6923 * s + 0.3077 * c;
    maps = {reshape(real(s), mapSize), reshape(real(c), mapSize), ...
        reshape(real(i), mapSize), reshape(real(e), mapSize)};
    limit = max(cellfun(@(map) max(abs(map(:))), maps));
    if limit <= eps; limit = 1; end
    for populationIndex = 1:4
        axisHandle = nexttile(layout);
        imagesc(axisHandle, maps{populationIndex});
        axis(axisHandle, 'image');
        axis(axisHandle, 'off');
        clim(axisHandle, [-limit limit]);
        colorbar(axisHandle, 'FontSize', 9);
        if modeIndex == 1
            title(axisHandle, populationLabels{populationIndex}, ...
                'FontSize', 12, 'FontWeight', 'bold');
        end
        if populationIndex == 1
            text(axisHandle, -0.08, 0.5, sprintf('mode %d: %.5g%+.5gi', modeIndex, ...
                real(modes.Lambda(modeIndex)), imag(modes.Lambda(modeIndex))), ...
                'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'middle', 'Rotation', 90, ...
                'Clipping', 'off', 'FontSize', 10);
        end
    end
end
title(layout, titleLines, 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'none');
l6ns_save_figure(figureHandle, outputRoot, stem);
close(figureHandle);
fprintf('Saved %s for angle %.2f, contrast %d, focus %s, w_I %s.\n', ...
    stem, angleValue, contrast, focus, local_display_number(wValue));
end

function [rows, stored] = local_mode_records(modes, angleValue, contrast, focus, ...
        stateLabel, wValue, gammaI, figureStem)
n = size(modes.Right, 1) / 3;
rows = cell(4, 16);
for modeIndex = 1:4
    vector = modes.Right(:, modeIndex);
    energy = [sum(abs(vector(1:n)).^2), ...
        sum(abs(vector(n + (1:n))).^2), ...
        sum(abs(vector(2*n + (1:n))).^2)];
    energy = energy / max(sum(energy), eps);
    rows(modeIndex, :) = {angleValue, contrast, string(focus), string(stateLabel), ...
        wValue, gammaI, modeIndex, real(modes.Lambda(modeIndex)), ...
        imag(modes.Lambda(modeIndex)), modes.ConditionNumber(modeIndex), ...
        modes.RightResidual(modeIndex), modes.LeftResidual(modeIndex), ...
        energy(1), energy(2), energy(3), string(figureStem)};
end
stored = struct('Lambda', modes.Lambda, 'Right', modes.Right, 'Left', modes.Left, ...
    'ConditionNumber', modes.ConditionNumber, ...
    'RightResidual', modes.RightResidual, 'LeftResidual', modes.LeftResidual, ...
    'WI', wValue, 'GammaI', gammaI, 'State', stateLabel, ...
    'FigureStem', figureStem);
end

function tag = local_angle_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end

function tag = local_number_tag(value, precision)
tag = sprintf(['%+.' num2str(precision) 'f'], value);
tag = strrep(tag, '+', 'p');
tag = strrep(tag, '-', 'm');
tag = strrep(tag, '.', 'p');
end

function textValue = local_display_number(value)
if isnan(value)
    textValue = 'component-only';
else
    textValue = sprintf('%+.9f', value);
end
end
