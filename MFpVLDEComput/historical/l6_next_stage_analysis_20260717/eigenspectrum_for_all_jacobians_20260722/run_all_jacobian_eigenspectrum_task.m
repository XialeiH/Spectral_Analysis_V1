function manifest = run_all_jacobian_eigenspectrum_task(taskIndex, ...
        l6SourceRoot, inhibitionSourceRoot, l6AtlasRoot, ...
        inhibitionAtlasRoot, outputRoot)
% Plot full-eigenvalue histograms and complex spectra for audited J families.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(l6SourceRoot)
    l6SourceRoot = getenv('ALLJ_L6_SOURCE_ROOT');
end
if nargin < 3 || isempty(inhibitionSourceRoot)
    inhibitionSourceRoot = getenv('ALLJ_INHIBITION_SOURCE_ROOT');
end
if nargin < 4 || isempty(l6AtlasRoot)
    l6AtlasRoot = getenv('ALLJ_L6_ATLAS_ROOT');
end
if nargin < 5 || isempty(inhibitionAtlasRoot)
    inhibitionAtlasRoot = getenv('ALLJ_INHIBITION_ATLAS_ROOT');
end
if nargin < 6 || isempty(outputRoot)
    outputRoot = getenv('ALLJ_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrast = 100;
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > 2 * numel(angles)
    error('AllJ:TaskIndex', 'Task index must be an integer from 1 to 8.');
end

if taskIndex <= numel(angles)
    family = 'L6';
    angleValue = angles(taskIndex);
    sourceRoot = l6SourceRoot;
    atlasRoot = l6AtlasRoot;
else
    family = 'Inhibition';
    angleValue = angles(taskIndex - numel(angles));
    sourceRoot = l6SourceRoot;
    atlasRoot = '';
end
local_require_folder(sourceRoot, 'source root');
if strcmp(family, 'L6')
    local_require_folder(atlasRoot, 'atlas root');
end
if isempty(outputRoot)
    error('AllJ:OutputRoot', 'ALLJ_OUTPUT_ROOT is required.');
end

if strcmp(family, 'L6')
    [j0, j1, eigenvalues0, eigenvalues1] = local_endpoint_jacobians( ...
        sourceRoot, angleValue, contrast);
    jBase = sparse(j1);
    jComponent = sparse(j0 - j1);
    [states, weights] = local_atlas_states(family, atlasRoot, angleValue, contrast);
else
    [j0, jBase, jComponent, eigenvalues0] = local_inhibition_jacobians( ...
        sourceRoot, angleValue, contrast);
    eigenvalues1 = eig(full(jBase), 'vector');
    [states, weights] = local_inhibition_states();
    [weights, boundaryAudit] = local_replace_boundary_weights( ...
        states, weights, jBase, jComponent);
end
componentReconstructionRelativeError = norm(j0 - (jBase + jComponent), 'fro') / ...
    max(norm(j0, 'fro'), eps);
if componentReconstructionRelativeError > 1e-12
    error('AllJ:Reconstruction', ...
        'Component reconstruction error is %.3e.', ...
        componentReconstructionRelativeError);
end
clear j0 j1
familyOutputRoot = fullfile(outputRoot, family);
if ~exist(familyOutputRoot, 'dir'); mkdir(familyOutputRoot); end

rows = cell(0, 17);
for stateIndex = 1:numel(weights)
    weight = weights(stateIndex);
    matrix = jBase + (1 - weight) * jComponent;
    if strcmp(family, 'L6')
        jacobianLabel = 'Total J(w)';
        definition = sprintf('J(w) = J_rest_6 + (1-w) J_6; w = %+.9f', weight);
        stem = sprintf('L6_Total_J_Angle_%sdeg_Contrast_%d_%s_w%s', ...
            local_angle_tag(angleValue), contrast, local_safe_name(states(stateIndex)), ...
            local_number_tag(weight, 9));
    else
        jacobianLabel = 'Total J(w_I)';
        definition = sprintf('J(wI) = J_rest_I + (1-wI) J_I; wI = %+.9f', weight);
        stem = sprintf('Inhibition_Total_J_Angle_%sdeg_Contrast_%d_%s_wI%s', ...
            local_angle_tag(angleValue), contrast, local_safe_name(states(stateIndex)), ...
            local_number_tag(weight, 9));
    end
    precomputedEigenvalues = [];
    if abs(weight) <= 1e-12
        precomputedEigenvalues = eigenvalues0;
    elseif abs(weight - 1) <= 1e-12
        precomputedEigenvalues = eigenvalues1;
    end
    [row, eigenvalues] = local_compute_and_plot(matrix, familyOutputRoot, stem, ...
        family, jacobianLabel, char(states(stateIndex)), definition, ...
        angleValue, contrast, weight, false, precomputedEigenvalues);
    rows(end + 1, :) = row; %#ok<AGROW>
    clear matrix eigenvalues
end

if strcmp(family, 'L6')
    jacobianLabel = 'Component J_6';
    stateLabel = 'L6 derivative component';
    definition = 'J_6 = J(w=0) - J(w=1)';
    stem = sprintf('L6_Component_J6_Angle_%sdeg_Contrast_%d', ...
        local_angle_tag(angleValue), contrast);
else
    jacobianLabel = 'Component J_I';
    stateLabel = 'Inhibitory-source Jacobian component';
    definition = 'J_I(:,I)=J_baseline(:,I); all S/C source columns are zero';
    stem = sprintf('Inhibition_Component_JI_Angle_%sdeg_Contrast_%d', ...
        local_angle_tag(angleValue), contrast);
end
[row, eigenvalues] = local_compute_and_plot(jComponent, familyOutputRoot, stem, ...
    family, jacobianLabel, stateLabel, definition, angleValue, contrast, ...
    NaN, true, []);
rows(end + 1, :) = row;
clear eigenvalues

if strcmp(family, 'L6')
    jacobianLabel = 'Rest Jacobian J_rest_6';
    stateLabel = 'L6-independent rest Jacobian';
    definition = 'J_rest_6 = J(w=1)';
    stem = sprintf('L6_Rest_Jrest6_Angle_%sdeg_Contrast_%d', ...
        local_angle_tag(angleValue), contrast);
else
    jacobianLabel = 'Rest Jacobian J_rest_I';
    stateLabel = 'Jacobian without inhibitory-source columns';
    definition = 'J_rest_I = J_baseline - J_I = J(wI=1)';
    stem = sprintf('Inhibition_Rest_JrestI_Angle_%sdeg_Contrast_%d', ...
        local_angle_tag(angleValue), contrast);
end
[row, eigenvalues] = local_compute_and_plot(jBase, familyOutputRoot, stem, ...
    family, jacobianLabel, stateLabel, definition, angleValue, contrast, ...
    NaN, true, eigenvalues1);
rows(end + 1, :) = row;
clear eigenvalues jBase jComponent

manifest = cell2table(rows, 'VariableNames', { ...
    'family','angle','contrast','jacobian','state','weight','eigenvalueCount', ...
    'binWidth','realAxisMin','realAxisMax','maxRealEigenvalue', ...
    'eigenRealMin','eigenRealMax','eigenImagMin','eigenImagMax', ...
    'histogramPdf','spectrumPdf'});
manifestFile = fullfile(outputRoot, sprintf( ...
    'manifest_%s_angle%s_contrast%d.tsv', lower(family), ...
    local_angle_tag(angleValue), contrast));
writetable(manifest, manifestFile, 'FileType', 'text', 'Delimiter', '\t');
if strcmp(family, 'Inhibition')
    boundaryFile = fullfile(outputRoot, sprintf( ...
        'boundary_%s_angle%s_contrast%d.tsv', lower(family), ...
        local_angle_tag(angleValue), contrast));
    writetable(boundaryAudit, boundaryFile, 'FileType', 'text', 'Delimiter', '\t');
end
fprintf('Saved %d %s Jacobian histogram/spectrum pairs for angle %.2f.\n', ...
    height(manifest), family, angleValue);
end

function [row, eigenvalues] = local_compute_and_plot(matrix, outputRoot, stem, ...
        family, jacobianLabel, stateLabel, definition, angleValue, contrast, ...
        weight, isComponent, precomputedEigenvalues)
if isempty(precomputedEigenvalues)
    fprintf('Computing full eigenspectrum: %s | angle %.2f | %s\n', ...
        jacobianLabel, angleValue, stateLabel);
    eigenvalues = eig(full(matrix), 'vector');
else
    fprintf('Using audited endpoint eigenspectrum: %s | angle %.2f | %s\n', ...
        jacobianLabel, angleValue, stateLabel);
    eigenvalues = precomputedEigenvalues(:);
end
if numel(eigenvalues) ~= size(matrix, 1) || any(~isfinite(eigenvalues))
    error('AllJ:Eigenvalues', 'Full eigenspectrum is incomplete or non-finite.');
end

binWidth = 0.05;
realValues = real(eigenvalues);
realMin = min(realValues);
realMax = max(realValues);
imagMin = min(imag(eigenvalues));
imagMax = max(imag(eigenvalues));
edgeMin = floor(realMin / binWidth) * binWidth;
edgeMax = ceil(realMax / binWidth) * binWidth;
if edgeMax <= edgeMin
    edgeMax = edgeMin + binWidth;
end
edges = edgeMin:binWidth:edgeMax;
if edges(end) < realMax
    edges(end + 1) = edges(end) + binWidth;
end
realLimits = [edges(1), edges(end)];
maxReal = max(realValues);

histogramPdf = fullfile(outputRoot, [stem, '_Histogram.pdf']);
spectrumPdf = fullfile(outputRoot, [stem, '_Eigenspectrum.pdf']);
local_plot_histogram(realValues, edges, realLimits, histogramPdf, ...
    jacobianLabel, stateLabel, definition, angleValue, contrast, maxReal, ...
    realMin, realMax);
local_plot_spectrum(eigenvalues, realLimits, spectrumPdf, ...
    jacobianLabel, stateLabel, definition, angleValue, contrast, maxReal, ...
    realMin, realMax, imagMin, imagMax);

if isComponent
    weightValue = NaN;
else
    weightValue = weight;
end
row = {family, angleValue, contrast, jacobianLabel, stateLabel, weightValue, ...
    numel(eigenvalues), binWidth, realLimits(1), realLimits(2), maxReal, ...
    realMin, realMax, imagMin, imagMax, ...
    histogramPdf, spectrumPdf};
end

function local_plot_histogram(realValues, edges, realLimits, outputFile, ...
        jacobianLabel, stateLabel, definition, angleValue, contrast, maxReal, ...
        realMin, realMax)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1400 950]);
ax = axes(fig);
ax.Position = [0.09 0.10 0.86 0.52];
histogram(ax, realValues, 'BinEdges', edges, ...
    'FaceColor', [0.0000 0.4470 0.7410], 'EdgeColor', [0.18 0.18 0.18], ...
    'LineWidth', 0.35);
hold(ax, 'on');
xline(ax, 0, ':', 'Color', [0.30 0.30 0.30], 'LineWidth', 1.1);
hasStabilityBoundary = realLimits(1) <= 1 && realLimits(2) >= 1;
if hasStabilityBoundary
    xline(ax, 1, '--', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.4);
end
xlim(ax, realLimits);
grid(ax, 'on');
box(ax, 'on');
ax.FontSize = 13;
xlabel(ax, 'Real part of eigenvalue, Re(lambda)', 'FontSize', 15);
ylabel(ax, 'Eigenvalue count', 'FontSize', 15);
title(ax, {sprintf('%s: eigenvalue histogram', jacobianLabel), ...
    sprintf('%s | angle %.2f deg | contrast %d', ...
        strrep(stateLabel, '_', ' '), angleValue, contrast), ...
    definition, ...
    sprintf('N=%d | eigenvalue Re range=[%.6g, %.6g]', ...
        numel(realValues), realMin, realMax), ...
    sprintf('bin width=0.05 | displayed Re axis=[%.2f, %.2f]', ...
        realLimits(1), realLimits(2)), ...
    sprintf('max Re(lambda)=%.6f', maxReal), ...
    sprintf('Reference lines: dotted gray Re(lambda)=0%s', ...
        local_stability_title(hasStabilityBoundary))}, ...
    'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');
exportgraphics(fig, outputFile, 'ContentType', 'vector');
close(fig);
end

function local_plot_spectrum(eigenvalues, realLimits, outputFile, ...
        jacobianLabel, stateLabel, definition, angleValue, contrast, maxReal, ...
        realMin, realMax, imagMin, imagMax)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1400 950]);
ax = axes(fig);
ax.Position = [0.09 0.10 0.86 0.52];
markerSize = 10;
scatter(ax, real(eigenvalues), imag(eigenvalues), markerSize, ...
    [0.0000 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.65, ...
    'MarkerEdgeAlpha', 0.25);
hold(ax, 'on');
xline(ax, 0, ':', 'Color', [0.30 0.30 0.30], 'LineWidth', 1.1);
hasStabilityBoundary = realLimits(1) <= 1 && realLimits(2) >= 1;
if hasStabilityBoundary
    xline(ax, 1, '--', 'Color', [0.85 0.15 0.15], 'LineWidth', 1.4);
end
yLimit = max(abs(imag(eigenvalues)));
if yLimit <= 1e-12
    yLimit = 0.05;
else
    yLimit = 1.04 * yLimit;
end
xlim(ax, realLimits);
ylim(ax, [-yLimit yLimit]);
grid(ax, 'on');
box(ax, 'on');
ax.FontSize = 13;
xlabel(ax, 'Real part, Re(lambda)', 'FontSize', 15);
ylabel(ax, 'Imaginary part, Im(lambda)', 'FontSize', 15);
title(ax, {sprintf('%s: complex eigenspectrum', jacobianLabel), ...
    sprintf('%s | angle %.2f deg | contrast %d', ...
        strrep(stateLabel, '_', ' '), angleValue, contrast), ...
    definition, ...
    sprintf('N=%d | eigenvalue Re range=[%.6g, %.6g]', ...
        numel(eigenvalues), realMin, realMax), ...
    sprintf('eigenvalue Im range=[%.6g, %.6g]', imagMin, imagMax), ...
    sprintf('displayed Re axis=[%.2f, %.2f] (matches histogram)', ...
        realLimits(1), realLimits(2)), ...
    sprintf('displayed Im axis=[%.2f, %.2f] | max Re(lambda)=%.6f', ...
        -yLimit, yLimit, maxReal), ...
    sprintf('Reference lines: dotted gray Re(lambda)=0%s', ...
        local_stability_title(hasStabilityBoundary))}, ...
    'Interpreter', 'none', 'FontSize', 11, 'FontWeight', 'bold');
exportgraphics(fig, outputFile, 'ContentType', 'vector');
close(fig);
end

function suffix = local_stability_title(hasStabilityBoundary)
if hasStabilityBoundary
    suffix = ' | dashed red: Re(lambda)=1';
else
    suffix = '';
end
end

function [j0, j1, eigenvalues0, eigenvalues1] = local_endpoint_jacobians( ...
        sourceRoot, angleValue, contrast)
file0 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW0p00_contr%d_angle_%.2f.mat', ...
    contrast, angleValue));
file1 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW1p00_contr%d_angle_%.2f.mat', ...
    contrast, angleValue));
local_require_file(file0, 'dynamic endpoint');
local_require_file(file1, 'frozen endpoint');
loaded0 = load(file0, 'Section4');
loaded1 = load(file1, 'Section4');
j0 = loaded0.Section4.A;
j1 = loaded1.Section4.A;
eigenvalues0 = loaded0.Section4.EigenValues(:);
eigenvalues1 = loaded1.Section4.EigenValues(:);
if ~isequal(size(j0), size(j1)) || size(j0, 1) ~= size(j0, 2)
    error('AllJ:JacobianSize', 'Endpoint Jacobians have incompatible dimensions.');
end
if numel(eigenvalues0) ~= size(j0, 1) || numel(eigenvalues1) ~= size(j1, 1)
    error('AllJ:EndpointEigenvalues', ...
        'Audited endpoint eigenspectra do not match the Jacobian dimensions.');
end
end

function [j0, jRestI, jI, eigenvalues0] = local_inhibition_jacobians( ...
        sourceRoot, angleValue, contrast)
file0 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW0p00_contr%d_angle_%.2f.mat', ...
    contrast, angleValue));
local_require_file(file0, 'baseline dynamic-L6 Jacobian');
loaded = load(file0, 'Section4');
j0 = sparse(loaded.Section4.A);
eigenvalues0 = loaded.Section4.EigenValues(:);
if size(j0, 1) ~= size(j0, 2) || mod(size(j0, 1), 3) ~= 0
    error('AllJ:InhibitionDimension', ...
        'Expected a square three-population baseline Jacobian.');
end
n = size(j0, 1) / 3;
iColumns = 2 * n + (1:n);
jI = sparse(size(j0, 1), size(j0, 2));
jI(:, iColumns) = j0(:, iColumns);
jRestI = sparse(j0 - jI);
end

function [states, weights] = local_inhibition_states()
states = ["Baseline_Dynamic_Inhibition_wI0"; ...
    "Frozen_Inhibition_wI1"; ...
    "Partial_Freeze_wI0p5"; ...
    "Highly_Stable_Stronger_Inhibition"; ...
    "Boundary_Just_Stable"; ...
    "Boundary_Just_Unstable"];
weights = [0; 1; 0.5; -1; NaN; NaN];
end

function [weights, audit] = local_replace_boundary_weights( ...
        states, weights, jBase, jComponent)
stableWeight = 0;
unstableWeight = 0.25;
while local_max_real(jBase + (1 - unstableWeight) * jComponent) <= 1
    unstableWeight = 2 * unstableWeight;
    if unstableWeight > 8
        error('AllJ:InhibitionBoundary', ...
            'Could not bracket the inhibition stability boundary.');
    end
end
for iteration = 1:30
    midpoint = 0.5 * (stableWeight + unstableWeight);
    if local_max_real(jBase + (1 - midpoint) * jComponent) > 1
        unstableWeight = midpoint;
    else
        stableWeight = midpoint;
    end
end
stableIndex = find(states == "Boundary_Just_Stable", 1);
unstableIndex = find(states == "Boundary_Just_Unstable", 1);
weights(stableIndex) = stableWeight;
weights(unstableIndex) = unstableWeight;
stableMaxReal = local_max_real(jBase + (1 - stableWeight) * jComponent);
unstableMaxReal = local_max_real(jBase + (1 - unstableWeight) * jComponent);
audit = table("Inhibition", stableWeight, unstableWeight, stableMaxReal, ...
    unstableMaxReal, abs(unstableWeight - stableWeight), ...
    'VariableNames', {'family','stableWeight','unstableWeight', ...
    'stableMaxReal','unstableMaxReal','bracketWidth'});
end

function value = local_max_real(matrix)
options = struct('tol', 1e-10, 'maxit', 3000, ...
    'p', min(60, size(matrix, 1)), 'isreal', true, ...
    'v0', ones(size(matrix, 1), 1) / sqrt(size(matrix, 1)));
try
    [~, diagonal, flag] = eigs(matrix, 12, 'largestreal', options);
    if flag ~= 0
        error('AllJ:EigsConvergence', 'eigs returned flag %d.', flag);
    end
    value = max(real(diag(diagonal)));
catch exception
    warning('AllJ:EigsFallback', ...
        'Rightmost eigs failed (%s); using full eig.', exception.message);
    value = max(real(eig(full(matrix), 'vector')));
end
end

function [states, weights] = local_atlas_states(family, atlasRoot, angleValue, contrast)
if strcmp(family, 'L6')
    summaryFile = fullfile(atlasRoot, 'results', 'data', sprintf( ...
        'top4_mode_summary_angle%s_contrast%d.tsv', ...
        local_angle_tag(angleValue), contrast));
    focus = 'Full_Jw';
    weightColumn = 'w';
else
    summaryFile = fullfile(atlasRoot, 'results', 'data', sprintf( ...
        'top4_inhibition_mode_summary_angle%s_contrast%d.tsv', ...
        local_angle_tag(angleValue), contrast));
    focus = 'Full_JwI';
    weightColumn = 'wI';
end
local_require_file(summaryFile, 'mode-atlas summary');
tableIn = readtable(summaryFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'VariableNamingRule', 'preserve');
rows = tableIn(string(tableIn.jacobianFocus) == focus, :);
[~, firstRows] = unique(string(rows.state), 'stable');
rows = rows(firstRows, :);
states = string(rows.state);
weights = rows.(weightColumn);
if numel(weights) ~= 6 || any(~isfinite(weights))
    error('AllJ:AtlasStates', ...
        'Expected six finite total-J weights in %s.', summaryFile);
end
end

function local_require_folder(pathValue, description)
if isempty(pathValue) || ~isfolder(pathValue)
    error('AllJ:MissingFolder', 'Missing %s: %s', description, pathValue);
end
end

function local_require_file(pathValue, description)
if ~isfile(pathValue)
    error('AllJ:MissingFile', 'Missing %s: %s', description, pathValue);
end
end

function tag = local_angle_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end

function tag = local_number_tag(value, digits)
tag = sprintf(['%+.', num2str(digits), 'f'], value);
tag = strrep(strrep(strrep(tag, '+', 'p'), '-', 'm'), '.', 'p');
end

function name = local_safe_name(value)
name = regexprep(char(value), '[^A-Za-z0-9]+', '_');
name = regexprep(name, '^_+|_+$', '');
end
