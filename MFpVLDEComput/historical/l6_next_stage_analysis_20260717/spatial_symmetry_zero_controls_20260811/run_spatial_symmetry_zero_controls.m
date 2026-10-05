function result = run_spatial_symmetry_zero_controls(sourceFile, canonicalFile, ...
        outputRoot, randomSeed, makeFigures)
% Test whether near-zero eigenvalues require spatial symmetry beyond block means.

if nargin < 1 || isempty(sourceFile)
    sourceFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'results_two_pathway_final', 'two_pathway_plan', ...
        'two_pathway_result.mat');
end
if nargin < 2 || isempty(canonicalFile)
    canonicalFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'near_zero_clustering_20260722', 'results_canonical', ...
        'near_zero_canonical_result.mat');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = fullfile(fileparts(mfilename('fullpath')), 'results');
end
if nargin < 4 || isempty(randomSeed); randomSeed = 8112026; end
if nargin < 5 || isempty(makeFigures); makeFigures = true; end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

loaded = load(sourceFile, 'result');
pathway = loaded.result.Pathway;
jBaseline = sparse(pathway.JBaseline);
n = pathway.PopulationSize;
mapSide = round(sqrt(n));
if mapSide^2 ~= n
    error('SpatialSymmetry:Map', 'Population size %d is not a square map.', n);
end
mapSize = [mapSide mapSide];
dimension = size(jBaseline, 1);
if dimension ~= 3*n
    error('SpatialSymmetry:Dimension', 'Expected a three-population Jacobian.');
end

canonical = load(canonicalFile, 'result');
baselineEigenvalues = canonical.result.Spectra.J_baseline(:);
valueShuffleEigenvalues = ...
    canonical.result.RandomSpectra.random_pattern_J_baseline(:);

rng(randomSeed, 'twister');
translationOnlyOrder = local_shift_order(mapSize, 3, 7);
d2Order = local_d2_orbit_permutation(mapSize);
sharedRandomOrder = randperm(n);
independentOrders = {randperm(n), randperm(n), randperm(n)};

sharedTranslationColumns = local_population_columns(translationOnlyOrder, n);
sharedD2Columns = local_population_columns(d2Order, n);
sharedRandomColumns = local_population_columns(sharedRandomOrder, n);
independentColumns = [independentOrders{1}, ...
    n + independentOrders{2}, 2*n + independentOrders{3}];

controlNames = { ...
    'baseline'; ...
    'preserve_translation_break_reflection_C2'; ...
    'preserve_reflection_C2_break_translation'; ...
    'break_translation_and_reflection_shared_SC'; ...
    'break_both_and_SC_alignment'; ...
    'shuffle_values_preserve_zero_mask'; ...
    'block_mean_only'};

controls = cell(numel(controlNames), 1);
controls{1} = jBaseline;
controls{2} = jBaseline(:, sharedTranslationColumns);
controls{3} = jBaseline(:, sharedD2Columns);
controls{4} = jBaseline(:, sharedRandomColumns);
controls{5} = jBaseline(:, independentColumns);
rng(7722, 'twister');
controls{6} = local_shuffle_blocks(jBaseline, n);
controls{7} = [];

transformNames = {'translate_row10', 'translate_column20', ...
    'reflect_rows', 'reflect_columns', 'rotate_180', 'rotate_90'};
transformOrders = { ...
    local_shift_order(mapSize, 10, 0), ...
    local_shift_order(mapSize, 0, 20), ...
    local_transform_order(mapSize, 'reflect_rows'), ...
    local_transform_order(mapSize, 'reflect_columns'), ...
    local_transform_order(mapSize, 'rotate_180'), ...
    local_transform_order(mapSize, 'rotate_90')};

symmetryRows = cell(numel(controlNames)*numel(transformNames), 4);
invariantRows = cell(numel(controlNames), 8);
spectra = cell(numel(controlNames), 1);
spectrumRows = cell(numel(controlNames), 13);
baselineBlockMeans = local_block_means(jBaseline, n);
baselineRowSums = full(sum(jBaseline, 2));

for controlIndex = 1:numel(controlNames)
    name = controlNames{controlIndex};
    fprintf('Control %d/%d: %s\n', controlIndex, numel(controlNames), name);
    if strcmp(name, 'block_mean_only')
        blockMeans = baselineBlockMeans;
        reducedMeanMatrix = n * blockMeans;
        activeEigenvalues = eig(reducedMeanMatrix, 'vector');
        eigenvalues = [activeEigenvalues; zeros(dimension - numel(activeEigenvalues), 1)];
        symmetryErrors = zeros(numel(transformNames), 1);
        meanError = 0;
        rowSumError = NaN;
        entryMultisetPreserved = false;
        zeroMaskPreserved = false;
        scResidual = local_mean_sc_residual(blockMeans);
        relativeMatrixChange = NaN;
    else
        matrix = controls{controlIndex};
        symmetryErrors = local_symmetry_errors(matrix, transformOrders, n);
        meanError = max(abs(local_block_means(matrix, n) - baselineBlockMeans), [], 'all');
        rowSumError = norm(full(sum(matrix, 2)) - baselineRowSums) / ...
            max(norm(baselineRowSums), eps);
        entryMultisetPreserved = controlIndex <= 5;
        zeroMaskPreserved = controlIndex == 1 || controlIndex == 6;
        scResidual = local_sc_proportionality_residual(matrix, n);
        relativeMatrixChange = norm(matrix - jBaseline, 'fro') / ...
            max(norm(jBaseline, 'fro'), eps);

        if controlIndex == 1
            eigenvalues = baselineEigenvalues;
        elseif controlIndex == 6
            eigenvalues = valueShuffleEigenvalues;
        else
            eigenvalues = local_full_eigenvalues(matrix);
        end
    end
    spectra{controlIndex} = eigenvalues;

    for symmetryIndex = 1:numel(transformNames)
        row = (controlIndex - 1)*numel(transformNames) + symmetryIndex;
        symmetryRows(row, :) = {string(name), string(transformNames{symmetryIndex}), ...
            symmetryErrors(symmetryIndex), ...
            symmetryErrors(symmetryIndex) <= 1e-10};
    end
    invariantRows(controlIndex, :) = {string(name), meanError, rowSumError, ...
        entryMultisetPreserved, zeroMaskPreserved, scResidual, ...
        relativeMatrixChange, numel(eigenvalues)};
    spectrumRows(controlIndex, :) = local_spectrum_row(name, eigenvalues);
end

symmetryAudit = cell2table(symmetryRows, 'VariableNames', ...
    {'control', 'symmetry', 'relativeCommutatorError', 'preservedAt1e10'});
invariantAudit = cell2table(invariantRows, 'VariableNames', ...
    {'control', 'maximumBlockMeanAbsError', 'relativeRowSumError', ...
    'blockEntryMultisetPreserved', 'zeroMaskPreserved', ...
    'SCProportionalityResidual', 'relativeMatrixChange', 'eigenvalueCount'});
spectrumSummary = cell2table(spectrumRows, 'VariableNames', ...
    {'control', 'maximumRealEigenvalue', 'minimumRealEigenvalue', ...
    'spectralRadius', 'minimumAbsEigenvalue', 'countAbsEigLE1e10', ...
    'countAbsEigLE1e4', 'countAbsEigLE1e2', 'countAbsEigLE2p5e2', ...
    'countAbsEigLE5e2', 'countAbsRealLE1e2', ...
    'countAbsRealLE2p5e2', 'countAbsRealLE5e2'});

local_write(symmetryAudit, outputRoot, 'symmetry_audit.tsv');
local_write(invariantAudit, outputRoot, 'control_invariants.tsv');
local_write(spectrumSummary, outputRoot, 'zero_count_summary.tsv');
if makeFigures
    local_plot_symmetry(symmetryAudit, controlNames, transformNames, outputRoot);
    local_plot_zero_counts(spectrumSummary, outputRoot);
    local_plot_histograms(controlNames, spectra, outputRoot);
end

result = struct('ControlNames', {controlNames}, 'SymmetryAudit', symmetryAudit, ...
    'InvariantAudit', invariantAudit, 'SpectrumSummary', spectrumSummary, ...
    'Spectra', {spectra}, 'SourceFile', sourceFile, ...
    'CanonicalFile', canonicalFile, 'MapSize', mapSize, 'RandomSeed', randomSeed);
save(fullfile(outputRoot, 'spatial_symmetry_zero_controls.mat'), ...
    'result', '-v7.3');
fprintf('Saved spatial-symmetry zero controls to %s\n', outputRoot);
end

function columns = local_population_columns(order, n)
columns = [order, n + order, 2*n + order];
end

function order = local_shift_order(mapSize, rowShift, columnShift)
indexMap = reshape(1:prod(mapSize), mapSize);
order = reshape(circshift(indexMap, [rowShift columnShift]), 1, []);
end

function order = local_transform_order(mapSize, transformName)
indexMap = reshape(1:prod(mapSize), mapSize);
switch transformName
    case 'reflect_rows'
        transformed = flipud(indexMap);
    case 'reflect_columns'
        transformed = fliplr(indexMap);
    case 'rotate_180'
        transformed = rot90(indexMap, 2);
    case 'rotate_90'
        transformed = rot90(indexMap, 1);
    otherwise
        error('Unknown transform %s.', transformName);
end
order = reshape(transformed, 1, []);
end

function order = local_d2_orbit_permutation(mapSize)
if any(mod(mapSize, 2) ~= 0)
    error('D2 orbit permutation expects even map dimensions.');
end
indexMap = reshape(1:prod(mapSize), mapSize);
halfRows = mapSize(1)/2;
halfColumns = mapSize(2)/2;
representatives = zeros(halfRows*halfColumns, 2);
entry = 0;
for column = 1:halfColumns
    for row = 1:halfRows
        entry = entry + 1;
        representatives(entry, :) = [row column];
    end
end
orbitPermutation = randperm(size(representatives, 1));
order = zeros(1, prod(mapSize));
for orbit = 1:size(representatives, 1)
    outputRepresentative = representatives(orbit, :);
    inputRepresentative = representatives(orbitPermutation(orbit), :);
    outputOrbit = local_d2_orbit(indexMap, outputRepresentative(1), ...
        outputRepresentative(2));
    inputOrbit = local_d2_orbit(indexMap, inputRepresentative(1), ...
        inputRepresentative(2));
    order(outputOrbit) = inputOrbit;
end
if ~isequal(sort(order), 1:prod(mapSize))
    error('D2 orbit construction did not produce a permutation.');
end
end

function orbit = local_d2_orbit(indexMap, row, column)
rowMirror = size(indexMap, 1) + 1 - row;
columnMirror = size(indexMap, 2) + 1 - column;
orbit = [indexMap(row, column), indexMap(rowMirror, column), ...
    indexMap(row, columnMirror), indexMap(rowMirror, columnMirror)];
end

function errors = local_symmetry_errors(matrix, transformOrders, n)
errors = zeros(numel(transformOrders), 1);
denominator = max(norm(matrix, 'fro'), eps);
for transform = 1:numel(transformOrders)
    populationOrder = local_population_columns(transformOrders{transform}, n);
    errors(transform) = norm(matrix(populationOrder, populationOrder) - matrix, 'fro') / ...
        denominator;
end
end

function means = local_block_means(matrix, n)
means = zeros(3, 3);
for target = 1:3
    rows = (target - 1)*n + (1:n);
    for source = 1:3
        columns = (source - 1)*n + (1:n);
        means(target, source) = full(sum(matrix(rows, columns), 'all')) / n^2;
    end
end
end

function residual = local_sc_proportionality_residual(matrix, n)
sColumns = matrix(:, 1:n);
cColumns = matrix(:, n + (1:n));
denominator = full(sum(sum(sColumns .* sColumns)));
if denominator <= eps
    residual = NaN;
    return
end
q = full(sum(sum(sColumns .* cColumns))) / denominator;
residual = norm(cColumns - q*sColumns, 'fro') / max(norm(cColumns, 'fro'), eps);
end

function residual = local_mean_sc_residual(blockMeans)
sColumn = blockMeans(:, 1);
cColumn = blockMeans(:, 2);
q = (sColumn' * cColumn) / max(sColumn' * sColumn, eps);
residual = norm(cColumn - q*sColumn) / max(norm(cColumn), eps);
end

function shuffled = local_shuffle_blocks(matrix, n)
shuffled = sparse(size(matrix, 1), size(matrix, 2));
for target = 1:3
    rows = (target - 1)*n + (1:n);
    for source = 1:3
        columns = (source - 1)*n + (1:n);
        block = matrix(rows, columns);
        [rowIndex, columnIndex, values] = find(block);
        if ~isempty(values)
            values = values(randperm(numel(values)));
            shuffled(rows, columns) = sparse(rowIndex, columnIndex, values, n, n);
        end
    end
end
end

function eigenvalues = local_full_eigenvalues(matrix)
denseMatrix = full(matrix);
eigenvalues = eig(denseMatrix, 'vector');
end

function row = local_spectrum_row(name, eigenvalues)
row = {string(name), max(real(eigenvalues)), min(real(eigenvalues)), ...
    max(abs(eigenvalues)), min(abs(eigenvalues)), ...
    sum(abs(eigenvalues) <= 1e-10), sum(abs(eigenvalues) <= 1e-4), ...
    sum(abs(eigenvalues) <= 1e-2), sum(abs(eigenvalues) <= 2.5e-2), ...
    sum(abs(eigenvalues) <= 5e-2), ...
    sum(abs(real(eigenvalues)) <= 1e-2), ...
    sum(abs(real(eigenvalues)) <= 2.5e-2), ...
    sum(abs(real(eigenvalues)) <= 5e-2)};
end

function local_plot_symmetry(audit, controlNames, transformNames, outputRoot)
values = nan(numel(controlNames), numel(transformNames));
for control = 1:numel(controlNames)
    for symmetry = 1:numel(transformNames)
        mask = audit.control == string(controlNames{control}) & ...
            audit.symmetry == string(transformNames{symmetry});
        values(control, symmetry) = audit.relativeCommutatorError(mask);
    end
end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 850]);
imagesc(log10(max(values, 1e-16)));
axis tight; colorbar; colormap(parula);
xticks(1:numel(transformNames)); xticklabels(local_pretty(transformNames));
yticks(1:numel(controlNames)); yticklabels(local_pretty(controlNames));
xtickangle(25);
xlabel('Tested spatial symmetry');
ylabel('Connectivity control');
title({'Spatial-symmetry audit of matched-mean Jacobian controls', ...
    'color = log_{10}(||PJP^{-1}-J||_F / ||J||_F)'});
exportgraphics(fig, fullfile(outputRoot, '01_spatial_symmetry_audit.pdf'), ...
    'ContentType', 'vector');
close(fig);
end

function local_plot_zero_counts(summary, outputRoot)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 850]);
counts = [summary.countAbsEigLE1e10, summary.countAbsEigLE1e2, ...
    summary.countAbsEigLE2p5e2, summary.countAbsEigLE5e2];
bar(counts, 'grouped');
grid on; box on;
xticks(1:height(summary)); xticklabels(local_pretty(cellstr(summary.control)));
xtickangle(24);
ylabel('Eigenvalue count');
title('Exact and near-zero eigenvalue counts after targeted spatial permutations');
legend({'|\lambda|\leq10^{-10}', '|\lambda|\leq0.01', ...
    '|\lambda|\leq0.025', '|\lambda|\leq0.05'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
exportgraphics(fig, fullfile(outputRoot, '02_exact_and_near_zero_counts.pdf'), ...
    'ContentType', 'vector');
close(fig);
end

function local_plot_histograms(controlNames, spectra, outputRoot)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1900 1050]);
layout = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for control = 1:numel(controlNames)
    values = real(spectra{control});
    edges = local_edges(values, 0.05);
    ax = nexttile(layout);
    histogram(ax, values, edges, 'FaceColor', [0.20 0.52 0.72], ...
        'EdgeColor', 'none');
    xline(ax, 0, 'k:');
    xlim(ax, [edges(1) edges(end)]);
    grid(ax, 'on');
    xlabel(ax, 'Re(\lambda)'); ylabel(ax, 'Count');
    title(ax, sprintf('%s\nN_{0.05}=%d, exact=%d', ...
        local_pretty(controlNames(control)), ...
        sum(abs(spectra{control}) <= 0.05), ...
        sum(abs(spectra{control}) <= 1e-10)), ...
        'FontWeight', 'normal', 'FontSize', 10);
end
title(layout, 'Full-Jacobian spectra under targeted spatial-structure controls');
exportgraphics(fig, fullfile(outputRoot, '03_control_eigenvalue_histograms.pdf'), ...
    'ContentType', 'vector');
close(fig);
end

function labels = local_pretty(labels)
labels = string(labels);
for index = 1:numel(labels)
    switch labels(index)
        case "baseline"
            labels(index) = "baseline";
        case "preserve_translation_break_reflection_C2"
            labels(index) = "translation kept; reflection/C2 broken";
        case "preserve_reflection_C2_break_translation"
            labels(index) = "reflection/C2 kept; translation broken";
        case "break_translation_and_reflection_shared_SC"
            labels(index) = "translation/reflection broken; S-C kept";
        case "break_both_and_SC_alignment"
            labels(index) = "translation/reflection/S-C broken";
        case "shuffle_values_preserve_zero_mask"
            labels(index) = "value shuffle; zero mask kept";
        case "block_mean_only"
            labels(index) = "3x3 block means only";
        case "translate_row10"
            labels(index) = "row translation 10";
        case "translate_column20"
            labels(index) = "column translation 20";
        case "reflect_rows"
            labels(index) = "row reflection";
        case "reflect_columns"
            labels(index) = "column reflection";
        case "rotate_180"
            labels(index) = "180-degree rotation";
        case "rotate_90"
            labels(index) = "90-degree rotation";
        otherwise
            labels(index) = replace(labels(index), '_', ' ');
    end
end
end

function edges = local_edges(values, width)
low = floor(min(values)/width)*width - width;
high = ceil(max(values)/width)*width + width;
edges = low:width:high;
end

function local_write(tableValue, outputRoot, filename)
writetable(tableValue, fullfile(outputRoot, filename), ...
    'FileType', 'text', 'Delimiter', '\t');
end
