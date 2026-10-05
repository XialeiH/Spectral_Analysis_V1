function summary = run_matched_destabilization_multicondition_task(taskIndex, sourceRoot, outputRoot)
% Repeat the matched L6/inhibition destabilization comparison across conditions.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('MATCHED_TASK_INDEX'));
    if isnan(taskIndex); taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID')); end
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('MATCHED_SOURCE_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('MATCHED_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
[contrastGrid, angleGrid] = ndgrid(contrasts, angles);
conditionCount = numel(contrastGrid);
assert(isscalar(taskIndex) && taskIndex == round(taskIndex) && ...
    taskIndex >= 1 && taskIndex <= conditionCount, ...
    'Task index must be an integer from 1 to %d.', conditionCount);
assert(isfolder(sourceRoot), 'Source root does not exist: %s', sourceRoot);
assert(~isempty(outputRoot), 'Output root is required.');

angleValue = angleGrid(taskIndex);
contrastValue = contrastGrid(taskIndex);
figureRoot = fullfile(outputRoot, 'figures');
dataRoot = fullfile(outputRoot, 'data');
if ~exist(figureRoot, 'dir'); mkdir(figureRoot); end
if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end

[j0, j1, sourceFiles] = load_endpoint_jacobians( ...
    sourceRoot, angleValue, contrastValue);
dimension = size(j0, 1);
assert(size(j0, 2) == dimension && isequal(size(j0), size(j1)), ...
    'Endpoint Jacobians must be square and equal in size.');
populationSize = dimension / 3;
side = round(sqrt(populationSize));
assert(populationSize == round(populationSize) && side^2 == populationSize, ...
    'Expected three square population maps.');

j6 = sparse(j0 - j1);
iColumns = 2*populationSize + (1:populationSize);
jI = sparse(dimension, dimension);
jI(:, iColumns) = j0(:, iColumns);
jRest = sparse(j0 - j6 - jI);
pathway = struct('JBaseline', sparse(j0), 'J6', j6, 'JI', jI, ...
    'JRest', jRest, 'PopulationSize', populationSize);
clear j0 j1

reconstructionError = norm( ...
    pathway.JRest + pathway.J6 + pathway.JI - pathway.JBaseline, 'fro') / ...
    max(norm(pathway.JBaseline, 'fro'), eps);
j6ISourceFraction = norm(pathway.J6(:, iColumns), 'fro') / ...
    max(norm(pathway.J6, 'fro'), eps);
jINonISourceFraction = norm(pathway.JI(:, 1:2*populationSize), 'fro') / ...
    max(norm(pathway.JI, 'fro'), eps);
assert(reconstructionError < 1e-12, ...
    'Pathway reconstruction error is %.3e.', reconstructionError);
assert(j6ISourceFraction < 1e-12 && jINonISourceFraction < 1e-12, ...
    'J6 and JI source-column definitions are not disjoint.');

rng(7000 + taskIndex, 'twister');
[lambdaBaseline, xBaseline, yBaseline] = critical_pair(pathway.JBaseline, 12);
[xBaseline, yBaseline] = orient_pair(xBaseline, yBaseline);
s6 = real(yBaseline' * pathway.J6 * xBaseline);
sI = real(yBaseline' * pathway.JI * xBaseline);
cStab = -s6 / sI;

delta6 = 0.05;
gamma6 = 1 + delta6;
jL6 = pathway.JRest + gamma6*pathway.J6 + pathway.JI;
alphaTarget = max_real(jL6);
[gammaI, matchMethod] = match_inhibition_gain(pathway, alphaTarget, cStab, delta6);
jInhibition = pathway.JRest + pathway.J6 + gammaI*pathway.JI;
alphaI = max_real(jInhibition);

[lambdaL6, xL6, yL6] = track_critical_branch( ...
    pathway, 'L6', gamma6, xBaseline, yBaseline);
[lambdaI, xI, yI] = track_critical_branch( ...
    pathway, 'I', gammaI, xBaseline, yBaseline);

weights = population_balance_weights(xL6, populationSize);
rightOverlap = abs((weights.*xL6)' * (weights.*xI)) / ...
    max(norm(weights.*xL6) * norm(weights.*xI), eps);
leftOverlap = abs((yL6./weights)' * (yI./weights)) / ...
    max(norm(yL6./weights) * norm(yI./weights), eps);

[rightL6, leftL6] = leading_subspaces(jL6, 4);
[rightI, leftI] = leading_subspaces(jInhibition, 4);
rightCosines = svd(orth(rightL6)' * orth(rightI));
leftCosines = svd(orth(leftL6)' * orth(leftI));

trackedGapL6 = alphaTarget - real(lambdaL6);
trackedGapI = alphaI - real(lambdaI);
alphaMatchError = abs(alphaTarget - alphaI);
expectedSigns = s6 > 0 && sI < 0;
matchedExactly = alphaMatchError < 1e-7;

summary = table(taskIndex, angleValue, contrastValue, real(lambdaBaseline), ...
    imag(lambdaBaseline), s6, sI, cStab, gamma6, gammaI, alphaTarget, alphaI, ...
    alphaMatchError, real(lambdaL6), real(lambdaI), trackedGapL6, trackedGapI, ...
    rightOverlap, leftOverlap, min(rightCosines), min(leftCosines), ...
    reconstructionError, j6ISourceFraction, jINonISourceFraction, ...
    expectedSigns, matchedExactly, string(matchMethod), ...
    string(sourceFiles.Dynamic), string(sourceFiles.Frozen), ...
    'VariableNames', {'taskIndex','angle','contrast','baselineLambdaReal', ...
    'baselineLambdaImag','s6','sI','cStab','gamma6','gammaI', ...
    'alphaL6','alphaI','alphaMatchError','trackedLambdaL6','trackedLambdaI', ...
    'trackedGapL6','trackedGapI','rightModeOverlap','leftModeOverlap', ...
    'minimumRightSubspaceCosine','minimumLeftSubspaceCosine', ...
    'reconstructionError','J6ISourceFraction','JINonISourceFraction', ...
    'expectedSensitivitySigns','matchedExactly','matchMethod', ...
    'dynamicEndpointFile','frozenEndpointFile'});

tag = sprintf('Angle_%sdeg_Contrast_%d', angle_tag(angleValue), contrastValue);
figureFile = fullfile(figureRoot, ...
    ['Matched_Destabilization_L6_Increase_vs_Inhibition_Decrease_Modes_', tag, '.pdf']);
plot_mode_comparison(xL6, xI, yL6, yI, populationSize, side, ...
    angleValue, contrastValue, gamma6, gammaI, alphaTarget, alphaI, ...
    lambdaL6, lambdaI, rightOverlap, leftOverlap, figureFile);

tableFile = fullfile(dataRoot, ['matched_destabilization_', tag, '.tsv']);
writetable(summary, tableFile, 'FileType', 'text', 'Delimiter', '\t');
save(fullfile(dataRoot, ['matched_destabilization_', tag, '.mat']), ...
    'summary', 'rightCosines', 'leftCosines', '-v7.3');
fprintf(['Completed angle %.2f, contrast %d: alpha %.9f vs %.9f, ', ...
    'right overlap %.6f, left overlap %.6f.\n'], ...
    angleValue, contrastValue, alphaTarget, alphaI, rightOverlap, leftOverlap);
end

function [j0, j1, files] = load_endpoint_jacobians(sourceRoot, angleValue, contrastValue)
dynamicTag = 'L6eqWS0p00_C0p00_I0p00';
frozenTag = 'L6eqWS1p00_C1p00_I1p00';
files.Dynamic = endpoint_file(sourceRoot, dynamicTag, angleValue, contrastValue);
files.Frozen = endpoint_file(sourceRoot, frozenTag, angleValue, contrastValue);
dynamic = load(files.Dynamic, 'Section4');
frozen = load(files.Frozen, 'Section4');
j0 = dynamic.Section4.A;
j1 = frozen.Section4.A;
end

function file = endpoint_file(sourceRoot, tag, angleValue, contrastValue)
file = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_%s_contr%d_angle_%.2f.mat', ...
    tag, contrastValue, angleValue));
assert(isfile(file), 'Missing endpoint file: %s', file);
end

function [gammaI, method] = match_inhibition_gain(pathway, alphaTarget, cStab, delta6)
objective = @(value) max_real( ...
    pathway.JRest + pathway.J6 + value*pathway.JI) - alphaTarget;
linearGuess = min(max(1 - cStab*delta6, 0), 1);
lower = max(0, linearGuess - 0.2);
upper = min(1, linearGuess + 0.2);
fLower = objective(lower);
fUpper = objective(upper);
if fLower * fUpper <= 0
    gammaI = fzero(objective, [lower upper]);
    method = 'local bracket';
    return
end
if lower > 0; fZero = objective(0); else; fZero = fLower; end
if upper < 1; fOne = objective(1); else; fOne = fUpper; end
if fZero * fOne <= 0
    gammaI = fzero(objective, [0 1]);
    method = 'full [0,1] bracket';
else
    gammaI = linearGuess;
    method = 'linear sensitivity fallback';
end
end

function [lambda, x, y] = track_critical_branch(pathway, axisName, gammaEnd, x, y)
gammaGrid = linspace(1, gammaEnd, 21);
lambda = NaN;
for step = 2:numel(gammaGrid)
    gamma = gammaGrid(step);
    if strcmp(axisName, 'L6')
        matrix = pathway.JRest + gamma*pathway.J6 + pathway.JI;
    else
        matrix = pathway.JRest + pathway.J6 + gamma*pathway.JI;
    end
    [lambdaPool, rightPool, leftPool] = paired_eigen_pool(matrix, 12);
    rightCosine = abs(y' * rightPool) ./ max(norm(y)*vecnorm(rightPool), eps);
    leftCosine = abs(leftPool' * x).' ./ max(vecnorm(leftPool)*norm(x), eps);
    score = sqrt(rightCosine .* leftCosine);
    [~, selected] = max(score);
    lambda = lambdaPool(selected);
    xNext = rightPool(:, selected);
    yNext = leftPool(:, selected);
    phase = angle(x' * xNext);
    x = xNext * exp(-1i*phase);
    y = yNext * exp(-1i*phase);
end
end

function [lambda, right, left] = paired_eigen_pool(matrix, count)
dimension = size(matrix, 1);
opts = struct('tol', 1e-9, 'maxit', 2000, ...
    'p', min(dimension, max(70, 3*count)));
[right, dRight] = eigs(matrix, count, 'largestreal', opts);
lambda = diag(dRight);
[leftCandidates, dLeft] = eigs(matrix', count + 6, 'largestreal', opts);
leftLambda = diag(dLeft);
left = zeros(size(right));
for k = 1:count
    right(:,k) = right(:,k) / norm(right(:,k));
    [~, index] = min(abs(leftLambda - conj(lambda(k))));
    candidate = leftCandidates(:,index);
    beta = candidate' * right(:,k);
    left(:,k) = candidate / conj(beta);
end
end

function [lambda, x, y] = critical_pair(matrix, poolSize)
dimension = size(matrix, 1);
opts = struct('tol', 1e-10, 'maxit', 2400, ...
    'p', min(dimension, max(60, 3*poolSize)));
[right, dRight] = eigs(matrix, poolSize, 'largestreal', opts);
lambdaPool = diag(dRight);
[~, index] = max(real(lambdaPool));
lambda = lambdaPool(index);
x = right(:, index);
[leftPool, dLeft] = eigs(matrix', poolSize, 'largestreal', opts);
leftLambda = diag(dLeft);
[~, leftIndex] = min(abs(leftLambda - conj(lambda)));
y = leftPool(:, leftIndex);
x = x / norm(x);
beta = y' * x;
assert(abs(beta) > 1e-10, 'Critical left/right vectors are nearly orthogonal.');
y = y / conj(beta);
end

function [x, y] = orient_pair(x, y)
[~, anchor] = max(abs(x));
phase = angle(x(anchor));
x = x * exp(-1i*phase);
y = y * exp(-1i*phase);
end

function alpha = max_real(matrix)
opts = struct('tol', 1e-9, 'maxit', 1800, 'p', 50);
values = eigs(matrix, 4, 'largestreal', opts);
alpha = max(real(values));
end

function weights = population_balance_weights(vector, n)
scale = [norm(vector(1:n))/sqrt(n), norm(vector(n+(1:n)))/sqrt(n), ...
    norm(vector(2*n+(1:n)))/sqrt(n)];
scale = max(scale, 1e-12);
weights = [ones(n,1)/scale(1); ones(n,1)/scale(2); ones(n,1)/scale(3)];
end

function [right, left] = leading_subspaces(matrix, count)
opts = struct('tol', 1e-9, 'maxit', 1800, 'p', 60);
[right, d] = eigs(matrix, count, 'largestreal', opts);
lambda = diag(d);
[leftPool, dLeft] = eigs(matrix', count + 4, 'largestreal', opts);
lambdaLeft = diag(dLeft);
left = zeros(size(matrix,1), count);
for k = 1:count
    [~, index] = min(abs(lambdaLeft - conj(lambda(k))));
    left(:,k) = leftPool(:,index);
end
end

function plot_mode_comparison(xL6, xI, yL6, yI, n, side, angleValue, ...
        contrastValue, gamma6, gammaI, alphaL6, alphaI, lambdaL6, lambdaI, ...
        rightOverlap, leftOverlap, file)
labels = {'S', 'C', 'I'};
indices = {1:n, n+(1:n), 2*n+(1:n)};
groups = {split_maps(real(xL6), indices, side), ...
    split_maps(real(xI), indices, side), ...
    split_maps(real(yL6), indices, side), ...
    split_maps(real(yI), indices, side)};
prefix = {sprintf('L6 increased (gamma_6=%.3f), right x', gamma6), ...
    sprintf('inhibition reduced (gamma_I=%.4f), right x', gammaI), ...
    'L6 increased, left y', 'inhibition reduced, left y'};

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1600 1600]);
t = tiledlayout(fig, 4, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
for row = 1:4
    for population = 1:3
        ax = nexttile(t);
        map = groups{row}{population};
        imagesc(ax, map);
        axis(ax, 'image');
        limit = max(max(abs(map(:))), eps);
        colormap(ax, parula(256));
        clim(ax, [-limit limit]);
        colorbar(ax, 'FontSize', 9);
        if row == 1
            title(ax, labels{population}, 'FontWeight', 'bold', 'FontSize', 12);
        end
        if row == 4
            xlabel(ax, 'spatial column');
        end
        if population == 1
            ylabel(ax, {prefix{row}, 'spatial row'}, ...
                'FontWeight', 'bold', 'FontSize', 10);
        end
    end
end
title(t, {sprintf(['Matched destabilization | angle %.2f deg | contrast %d | ', ...
    'full-J max Re lambda %.6f vs %.6f'], ...
    angleValue, contrastValue, alphaL6, alphaI), ...
    sprintf(['tracked Re lambda %.6f vs %.6f | weighted right overlap %.3f | ', ...
    'left overlap %.3f'], real(lambdaL6), real(lambdaI), rightOverlap, leftOverlap)}, ...
    'FontSize', 15, 'FontWeight', 'bold');
set(fig, 'Renderer', 'painters');
exportgraphics(fig, file, 'ContentType', 'vector');
close(fig);
end

function maps = split_maps(vector, indices, side)
maps = cell(1, numel(indices));
for k = 1:numel(indices)
    maps{k} = reshape(vector(indices{k}), side, side);
end
end

function tag = angle_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end
