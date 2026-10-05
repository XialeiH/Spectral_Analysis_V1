function summary = run_inhibition_revision_mode_task(taskIndex, sourceRoot, outputRoot)
% Plot matched eigenmodes and physical one-tau singular modes for inhibition.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('INHR_SOURCE_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('INHR_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
stateCount = 6;
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > numel(angles) * stateCount
    error('InhibitionRevision:ModeTask', ...
        'Task index must be an integer from 1 to 24.');
end
if isempty(sourceRoot) || ~isfolder(sourceRoot)
    error('InhibitionRevision:SourceRoot', ...
        'INHR_SOURCE_ROOT must identify the celltype-control run root.');
end
if isempty(outputRoot)
    error('InhibitionRevision:OutputRoot', 'INHR_OUTPUT_ROOT is required.');
end

angleIndex = floor((taskIndex - 1) / stateCount) + 1;
stateIndex = mod(taskIndex - 1, stateCount) + 1;
angle = angles(angleIndex);
contrast = 100;
boundaryFile = fullfile(outputRoot, 'boundary_data', sprintf( ...
    'inhibition_boundary_angle%.2f_contrast%d.tsv', angle, contrast));
if ~isfile(boundaryFile)
    error('InhibitionRevision:MissingBoundary', ...
        'Missing exact stability bracket %s.', boundaryFile);
end
boundary = readtable(boundaryFile, 'FileType', 'text', 'Delimiter', '\t');
wValues = [-0.25 0 boundary.lowW(1) boundary.highW(1) 0.5 1];
stateLabels = {'Stronger_Dynamic_Inhibition','Baseline_Dynamic_Inhibition', ...
    'Boundary_Just_Stable','Boundary_Just_Unstable','Partial_Freeze', ...
    'Frozen_Inhibition'};
wI = wValues(stateIndex);
stateLabel = stateLabels{stateIndex};

[jacobian0, jacobian1] = local_endpoint_jacobians(sourceRoot, angle, contrast);
jacobian = jacobian0 + wI * (jacobian1 - jacobian0);
clear jacobian0 jacobian1

cfg = l6ns_config();
cfg.EigsTolerance = 1e-10;
cfg.EigsMaxIterations = 2200;
cfg.EigsSubspaceDimension = 80;
cfg.TransientPowerIterations = 60;
cfg.TransientRefinePeak = true;
cfg.RandomSeed = 1700 + taskIndex;
subspaceOverride = str2double(getenv('INHR_TRANSIENT_SUBSPACE'));
if isfinite(subspaceOverride) && subspaceOverride >= 4
    cfg.TransientSvdsSubspaceDimension = round(subspaceOverride);
end

modes = l6ns_eigenpairs(jacobian, 6, 'largestreal', cfg);
lambda = modes.Lambda(1);
rightMode = modes.Right(:,1);
leftMode = modes.Left(:,1);
[rightMode, leftMode] = local_phase_pair(rightMode, leftMode);

tauMs = cfg.TauMs;
transient = l6ns_continuous_transient(jacobian, tauMs, cfg, [], tauMs);
[singularInput, singularOutput] = local_sign_pair( ...
    transient.OptimalInput, transient.PeakOutput);

stateSize = size(jacobian,1) / 3;
mapWidth = sqrt(stateSize);
if stateSize ~= round(stateSize) || mapWidth ~= round(mapWidth)
    error('InhibitionRevision:MapSize', ...
        'The 3-population state does not form square maps.');
end
mapSize = [mapWidth mapWidth];
components = {local_components(rightMode, stateSize), ...
    local_components(leftMode, stateSize), ...
    local_components(singularInput, stateSize), ...
    local_components(singularOutput, stateSize)};
rowNames = {'Right eigenmode','Matched left eigenmode', ...
    'One-tau singular input','One-tau singular output'};
populationNames = {'S','C','I','E = 0.6923 S + 0.3077 C'};

figureRoot = fullfile(outputRoot, 'mode_figures');
dataRoot = fullfile(outputRoot, 'mode_data');
if ~exist(figureRoot, 'dir'); mkdir(figureRoot); end
if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end

figureHandle = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [40 40 2800 2250]);
layout = tiledlayout(figureHandle, 4, 4, ...
    'TileSpacing', 'loose', 'Padding', 'loose');
colormap(figureHandle, parula(256));
for row = 1:4
    for column = 1:4
        axisHandle = nexttile(layout);
        map = reshape(real(components{row}{column}), mapSize);
        imagesc(axisHandle, map);
        axis(axisHandle, 'image');
        set(axisHandle, 'YDir', 'normal', 'FontSize', 11, 'LineWidth', 0.8);
        limit = max(abs(map(:)));
        if limit <= eps; limit = 1; end
        clim(axisHandle, [-limit limit]);
        colorbar(axisHandle, 'FontSize', 10);
        xlabel(axisHandle, 'cortical x index', 'FontSize', 11);
        ylabel(axisHandle, 'cortical y index', 'FontSize', 11);
        title(axisHandle, sprintf('%s: %s', rowNames{row}, ...
            populationNames{column}), 'FontSize', 12, 'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end
end
stableText = 'stable';
if real(lambda) >= 1; stableText = 'unstable'; end
titleLine1 = sprintf(['Inhibition mechanism modes | angle %.2f deg | ' ...
    'contrast %d | %s'], angle, contrast, strrep(stateLabel, '_', ' '));
titleLine2 = sprintf(['w_I=%+.9f | gamma_I=1-w_I=%.9f | ' ...
    'max Re(lambda)=%.9f | margin=%.9f (%s) | ' ...
    'sigma_1[exp((J-I)t/tau)] at t=tau=%.2f ms: %.6g'], ...
    wI, 1-wI, real(lambda), 1-real(lambda), stableText, ...
    tauMs, transient.PeakGain);
title(layout, {titleLine1, titleLine2}, 'FontSize', 16, ...
    'FontWeight', 'bold', 'Interpreter', 'none');

wTag = local_number_tag(wI, 6);
stem = sprintf('Inhibition_Modes_Angle_%.2fdeg_Contrast_%d_%s_wI%s', ...
    angle, contrast, stateLabel, wTag);
pdfFile = fullfile(figureRoot, [stem '.pdf']);
exportgraphics(figureHandle, pdfFile, 'ContentType', 'vector');
close(figureHandle);

summary = table(angle, contrast, wI, 1-wI, string(stateLabel), ...
    real(lambda), imag(lambda), 1-real(lambda), transient.PeakGain, tauMs, ...
    modes.ConditionNumber(1), modes.RightResidual(1), modes.LeftResidual(1), ...
    transient.ForwardSingularResidual, transient.AdjointSingularResidual, ...
    string(pdfFile), 'VariableNames', {'angle','contrast','wI','gammaI','state', ...
    'lambdaReal','lambdaImag','stabilityMargin','singularGainAtTau','tauMs', ...
    'eigenConditionNumber','rightEigenResidual','leftEigenResidual', ...
    'forwardSingularResidual','adjointSingularResidual','pdfFile'});
writetable(summary, fullfile(dataRoot, [stem '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
metadata = struct( ...
    'JacobianDefinition', 'J_I(w_I)=J_I(0)+w_I*(J_I(1)-J_I(0))', ...
    'PropagatorDefinition', 'exp((J-I)*t/tau)', 'HorizonMs', tauMs, ...
    'W0Meaning', 'original dynamic inhibition', ...
    'W1Meaning', 'frozen inhibitory derivative with dynamic L6');
save(fullfile(dataRoot, [stem '.mat']), 'summary', 'metadata', ...
    'rightMode', 'leftMode', 'singularInput', 'singularOutput', '-v7.3');
fprintf('Saved inhibition mode atlas to %s.\n', pdfFile);
end

function components = local_components(vector, stateSize)
s = vector(1:stateSize);
c = vector(stateSize+1:2*stateSize);
i = vector(2*stateSize+1:3*stateSize);
e = 0.6923 * s + 0.3077 * c;
components = {s,c,i,e};
end

function [rightMode, leftMode] = local_phase_pair(rightMode, leftMode)
[~, index] = max(abs(rightMode));
phase = angle(rightMode(index));
factor = exp(-1i * phase);
rightMode = rightMode * factor;
leftMode = leftMode * factor;
end

function [inputMode, outputMode] = local_sign_pair(inputMode, outputMode)
[~, index] = max(abs(inputMode));
if real(inputMode(index)) < 0
    inputMode = -inputMode;
    outputMode = -outputMode;
end
end

function tag = local_number_tag(value, precision)
tag = sprintf(['%+.' num2str(precision) 'f'], value);
tag = strrep(tag, '+', 'p');
tag = strrep(tag, '-', 'm');
tag = strrep(tag, '.', 'p');
end

function [jacobian0, jacobian1] = local_endpoint_jacobians(sourceRoot, angle, contrast)
file0 = local_endpoint_file(sourceRoot, angle, contrast, 0);
file1 = local_endpoint_file(sourceRoot, angle, contrast, 1);
loaded0 = load(file0, 'Section4');
loaded1 = load(file1, 'Section4');
jacobian0 = loaded0.Section4.A;
jacobian1 = loaded1.Section4.A;
if ~isequal(size(jacobian0), size(jacobian1))
    error('InhibitionRevision:JacobianSize', ...
        'Endpoint Jacobians have incompatible dimensions.');
end
end

function file = local_endpoint_file(sourceRoot, angle, contrast, wI)
tag = sprintf('L6eqWS0p00_C0p00_I%.2f', wI);
tag = strrep(strrep(tag, '.', 'p'), '-', 'm');
file = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_%s_contr%d_angle_%.2f.mat', ...
    tag, contrast, angle));
if ~isfile(file)
    error('InhibitionRevision:MissingEndpoint', 'Missing endpoint file %s.', file);
end
end
