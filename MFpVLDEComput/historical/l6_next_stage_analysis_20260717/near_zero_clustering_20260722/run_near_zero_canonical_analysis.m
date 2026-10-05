function result = run_near_zero_canonical_analysis(sourceFile, outputRoot)
% Detailed near-zero-cluster controls at angle 0 deg, contrast 100.

if nargin < 1 || isempty(sourceFile)
    sourceFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'results_two_pathway_final', 'two_pathway_plan', ...
        'two_pathway_result.mat');
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = fullfile(fileparts(mfilename('fullpath')), 'results_canonical');
end
if ~isfile(sourceFile)
    error('NearZeroCanonical:Source', 'Source file does not exist: %s', sourceFile);
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

loaded = load(sourceFile, 'result');
pathway = loaded.result.Pathway;
jRest = sparse(pathway.JRest);
j6 = sparse(pathway.J6);
jI = sparse(pathway.JI);
jBaseline = sparse(pathway.JBaseline);
n = pathway.PopulationSize;
dimension = size(jBaseline, 1);
eColumns = 1:2*n;
iColumns = 2*n + (1:n);

reconstructionError = norm(jRest + j6 + jI - jBaseline, 'fro') / ...
    max(norm(jBaseline, 'fro'), eps);
if reconstructionError > 1e-13
    error('NearZeroCanonical:Reconstruction', ...
        'J_rest + J_6 + J_I failed to reconstruct the baseline.');
end

matrixDefinitions = { ...
    'J_baseline', jBaseline; ...
    'J_rest', jRest; ...
    'J_6', j6; ...
    'J_I', jI; ...
    'J_rest_plus_J6', jRest + j6; ...
    'J_rest_plus_JI', jRest + jI; ...
    'J6_plus_JI', j6 + jI};

spectra = struct();
summaryRows = cell(size(matrixDefinitions, 1), 23);
for matrixIndex = 1:size(matrixDefinitions, 1)
    name = matrixDefinitions{matrixIndex, 1};
    matrix = sparse(matrixDefinitions{matrixIndex, 2});
    fprintf('Canonical spectrum: %s\n', name);
    eigenvalues = local_block_eigenvalues(matrix, n);
    spectra.(name) = eigenvalues;
    summaryRows(matrixIndex, :) = local_spectrum_row(name, matrix, ...
        eigenvalues, n, 'model');
end
matrixSummary = cell2table(summaryRows, 'VariableNames', local_summary_names());

rankRows = cell(4, 11);
rankRows(1, :) = local_rank_row('J_baseline', jBaseline, jBaseline, dimension, 0);
rankRows(2, :) = local_rank_row('J_rest', jRest, jRest(:, eColumns), 2*n, n);
[q6, residual6] = local_sc_proportionality(j6, n);
rankRows(3, :) = local_rank_row('J_6', j6, j6(:, 1:n), n, 2*n);
rankRows(4, :) = local_rank_row('J_I', jI, jI(:, iColumns), n, 2*n);
rankSummary = cell2table(rankRows, 'VariableNames', { ...
    'matrix', 'dimension', 'nnz', 'density', 'structuralRank', ...
    'independentSourceColumns', 'guaranteedGeometricNullity', ...
    'smallestActiveSingularValue', 'largestSingularValue', ...
    'stableRank', 'activeRankTolerance'});

gammaValues = [-4 -3 -2 -1 0 0.5 1 2 3 4];
proportionalityRows = cell(numel(gammaValues) + 2, 9);
rowIndex = 0;
for gamma6 = gammaValues
    rowIndex = rowIndex + 1;
    matrix = jRest + gamma6 * j6;
    proportionalityRows(rowIndex, :) = local_proportionality_row( ...
        sprintf('J_rest_plus_%+.3g_J6', gamma6), matrix, gamma6, n);
end
rowIndex = rowIndex + 1;
proportionalityRows(rowIndex, :) = local_proportionality_row( ...
    'J_rest', jRest, NaN, n);
rowIndex = rowIndex + 1;
proportionalityRows(rowIndex, :) = local_proportionality_row( ...
    'J_6', j6, NaN, n);
proportionalitySummary = cell2table(proportionalityRows, 'VariableNames', { ...
    'matrix', 'gamma6', 'bestScalarQ', 'globalRelativeResidual', ...
    'differenceActionFrobeniusNorm', 'sourceFrobeniusNorm', ...
    'differenceToSourceNormRatio', 'minimumColumnQ', 'maximumColumnQ'});

fprintf('Canonical full eigendecomposition for near-zero mode geometry.\n');
[rightVectors, baselineEigenvalues] = eig(full(jBaseline), 'vector');
spectra.J_baseline = baselineEigenvalues;
[qBaseline, ~] = local_sc_proportionality(jRest + j6, n);
thresholds = [0.01 0.025 0.05];
modeRows = cell(2 * numel(thresholds), 11);
modeRow = 0;
for threshold = thresholds
    nearMask = abs(baselineEigenvalues) <= threshold;
    for classIndex = 1:2
        modeRow = modeRow + 1;
        if classIndex == 1
            mask = nearMask;
            className = 'near_zero';
        else
            mask = ~nearMask;
            className = 'outside';
        end
        modeRows(modeRow, :) = local_mode_geometry_row(className, threshold, ...
            baselineEigenvalues(mask), rightVectors(:, mask), qBaseline, n);
    end
end
modeGeometry = cell2table(modeRows, 'VariableNames', { ...
    'class', 'absEigenvalueThreshold', 'modeCount', 'medianAbsEigenvalue', ...
    'meanSEnergyFraction', 'meanCEnergyFraction', 'meanIEnergyFraction', ...
    'medianSEnergyFraction', 'medianCEnergyFraction', ...
    'medianIEnergyFraction', 'medianSCCancellationScore'});
clear rightVectors

[schurSummary, schurSpectra] = local_schur_audit(jBaseline, n, ...
    baselineEigenvalues);

controlDefinitions = local_ablation_matrices(jRest, j6, jI, n);
controlRows = cell(size(controlDefinitions, 1), 23);
controlSpectra = struct();
for controlIndex = 1:size(controlDefinitions, 1)
    name = controlDefinitions{controlIndex, 1};
    matrix = sparse(controlDefinitions{controlIndex, 2});
    fprintf('Ablation spectrum: %s\n', name);
    eigenvalues = local_block_eigenvalues(matrix, n);
    controlSpectra.(name) = eigenvalues;
    controlRows(controlIndex, :) = local_spectrum_row(name, matrix, ...
        eigenvalues, n, 'ablation');
end
ablationSummary = cell2table(controlRows, 'VariableNames', local_summary_names());

rng(7722, 'twister');
randomDefinitions = { ...
    'random_pattern_J_baseline', local_shuffle_blocks(jBaseline, n); ...
    'random_pattern_J_rest', local_shuffle_blocks(jRest, n); ...
    'random_pattern_J_6', local_shuffle_blocks(j6, n); ...
    'random_pattern_J_I', local_shuffle_blocks(jI, n)};
randomRows = cell(size(randomDefinitions, 1), 23);
randomSpectra = struct();
for randomIndex = 1:size(randomDefinitions, 1)
    name = randomDefinitions{randomIndex, 1};
    matrix = sparse(randomDefinitions{randomIndex, 2});
    fprintf('Random-pattern control: %s\n', name);
    eigenvalues = local_block_eigenvalues(matrix, n);
    randomSpectra.(name) = eigenvalues;
    randomRows(randomIndex, :) = local_spectrum_row(name, matrix, ...
        eigenvalues, n, 'random_pattern_control');
end
randomControlSummary = cell2table(randomRows, ...
    'VariableNames', local_summary_names());

sweepGamma = [-4 -2 -1 0 0.5 1 2 4];
sweepDefinitions = cell(0, 3);
for gamma6 = sweepGamma
    sweepDefinitions(end + 1, :) = {'gamma6_sweep', gamma6, 1}; %#ok<AGROW>
end
for gammaI = sweepGamma
    sweepDefinitions(end + 1, :) = {'gammaI_sweep', 1, gammaI}; %#ok<AGROW>
end
jointGamma = [-2 0 1 3];
for gamma6 = jointGamma
    for gammaI = jointGamma
        sweepDefinitions(end + 1, :) = {'joint_grid', gamma6, gammaI}; %#ok<AGROW>
    end
end
sweepRows = cell(size(sweepDefinitions, 1), 17);
sweepSpectra = cell(size(sweepDefinitions, 1), 1);
for sweepIndex = 1:size(sweepDefinitions, 1)
    sweepType = sweepDefinitions{sweepIndex, 1};
    gamma6 = sweepDefinitions{sweepIndex, 2};
    gammaI = sweepDefinitions{sweepIndex, 3};
    matrix = jRest + gamma6*j6 + gammaI*jI;
    fprintf('Weight spectrum: %s gamma6=%+.3g gammaI=%+.3g\n', ...
        sweepType, gamma6, gammaI);
    eigenvalues = local_block_eigenvalues(matrix, n);
    sweepSpectra{sweepIndex} = eigenvalues;
    sweepRows(sweepIndex, :) = {sweepType, gamma6, gammaI, 1-gamma6, ...
        1-gammaI, max(real(eigenvalues)), min(real(eigenvalues)), ...
        max(abs(eigenvalues)), min(abs(eigenvalues)), ...
        sum(abs(eigenvalues) <= 1e-10), sum(abs(eigenvalues) <= 1e-4), ...
        sum(abs(eigenvalues) <= 1e-2), sum(abs(eigenvalues) <= 2.5e-2), ...
        sum(abs(eigenvalues) <= 5e-2), ...
        sum(abs(real(eigenvalues)) <= 1e-2), ...
        sum(abs(real(eigenvalues)) <= 2.5e-2), ...
        sum(abs(real(eigenvalues)) <= 5e-2)};
end
weightSweepSummary = cell2table(sweepRows, 'VariableNames', { ...
    'sweepType', 'gamma6', 'gammaI', 'w6', 'wI', ...
    'maximumRealEigenvalue', 'minimumRealEigenvalue', 'spectralRadius', ...
    'minimumAbsEigenvalue', 'countAbsEigLE1e10', 'countAbsEigLE1e4', ...
    'countAbsEigLE1e2', 'countAbsEigLE2p5e2', 'countAbsEigLE5e2', ...
    'countAbsRealLE1e2', 'countAbsRealLE2p5e2', 'countAbsRealLE5e2'});

audit = table(reconstructionError, q6, residual6, ...
    norm(j6(:, iColumns), 'fro') / max(norm(j6, 'fro'), eps), ...
    norm(jI(:, eColumns), 'fro') / max(norm(jI, 'fro'), eps), ...
    'VariableNames', {'reconstructionRelativeError', 'J6SCScalarQ', ...
    'J6SCProportionalityResidual', 'J6ISourceFraction', ...
    'JINonISourceFraction'});

local_write(matrixSummary, outputRoot, 'canonical_matrix_summary.tsv');
local_write(rankSummary, outputRoot, 'canonical_rank_summary.tsv');
local_write(proportionalitySummary, outputRoot, ...
    'canonical_sc_proportionality.tsv');
local_write(modeGeometry, outputRoot, 'canonical_near_zero_mode_geometry.tsv');
local_write(schurSummary, outputRoot, 'canonical_schur_summary.tsv');
local_write(ablationSummary, outputRoot, 'canonical_ablation_summary.tsv');
local_write(randomControlSummary, outputRoot, ...
    'canonical_random_pattern_controls.tsv');
local_write(weightSweepSummary, outputRoot, 'canonical_weight_sweep.tsv');
local_write(audit, outputRoot, 'canonical_audit.tsv');

result = struct('Audit', audit, 'MatrixSummary', matrixSummary, ...
    'RankSummary', rankSummary, 'Proportionality', proportionalitySummary, ...
    'ModeGeometry', modeGeometry, 'SchurSummary', schurSummary, ...
    'Ablations', ablationSummary, 'RandomControls', randomControlSummary, ...
    'WeightSweep', weightSweepSummary, 'Spectra', spectra, ...
    'ControlSpectra', controlSpectra, 'RandomSpectra', randomSpectra, ...
    'SweepSpectra', {sweepSpectra}, 'SchurSpectra', schurSpectra, ...
    'SourceFile', sourceFile);
save(fullfile(outputRoot, 'near_zero_canonical_result.mat'), 'result', '-v7.3');
fprintf('Canonical near-zero audit saved to %s.\n', outputRoot);
end

function row = local_rank_row(name, matrix, activeMatrix, independentColumns, ...
        guaranteedNullity)
largestSingular = svds(matrix, 1, 'largest');
smallestActive = svds(activeMatrix, 1, 'smallest', ...
    'Tolerance', 1e-9, 'MaxIterations', 3000);
rankTolerance = max(size(activeMatrix)) * eps(largestSingular);
row = {name, size(matrix, 1), nnz(matrix), nnz(matrix)/numel(matrix), ...
    sprank(matrix), independentColumns, guaranteedNullity, smallestActive, ...
    largestSingular, norm(matrix, 'fro')^2 / largestSingular^2, rankTolerance};
end

function [q, residual] = local_sc_proportionality(matrix, n)
sColumns = matrix(:, 1:n);
cColumns = matrix(:, n + (1:n));
denominator = full(sum(sum(sColumns .* sColumns)));
q = full(sum(sum(sColumns .* cColumns))) / denominator;
residual = norm(cColumns - q*sColumns, 'fro') / ...
    max(norm(cColumns, 'fro'), eps);
end

function row = local_proportionality_row(name, matrix, gamma6, n)
sColumns = matrix(:, 1:n);
cColumns = matrix(:, n + (1:n));
[q, residual] = local_sc_proportionality(matrix, n);
difference = cColumns - q*sColumns;
columnQ = nan(n, 1);
for column = 1:n
    source = sColumns(:, column);
    target = cColumns(:, column);
    denominator = full(source' * source);
    if denominator > eps
        columnQ(column) = full(source' * target) / denominator;
    end
end
valid = columnQ(isfinite(columnQ));
row = {name, gamma6, q, residual, norm(difference, 'fro'), ...
    norm([sColumns cColumns], 'fro'), ...
    norm(difference, 'fro') / max(norm([sColumns cColumns], 'fro'), eps), ...
    min(valid), max(valid)};
end

function row = local_mode_geometry_row(className, threshold, eigenvalues, ...
        vectors, q, n)
count = numel(eigenvalues);
if count == 0
    row = {className, threshold, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN};
    return
end
sEnergy = zeros(count, 1);
cEnergy = zeros(count, 1);
iEnergy = zeros(count, 1);
cancellation = zeros(count, 1);
for mode = 1:count
    s = vectors(1:n, mode);
    c = vectors(n + (1:n), mode);
    i = vectors(2*n + (1:n), mode);
    energy = [sum(abs(s).^2), sum(abs(c).^2), sum(abs(i).^2)];
    energy = energy / max(sum(energy), eps);
    sEnergy(mode) = energy(1);
    cEnergy(mode) = energy(2);
    iEnergy(mode) = energy(3);
    cancellation(mode) = norm(s + q*c) / ...
        max(sqrt(norm(s)^2 + q^2*norm(c)^2), eps);
end
row = {className, threshold, count, median(abs(eigenvalues)), ...
    mean(sEnergy), mean(cEnergy), mean(iEnergy), median(sEnergy), ...
    median(cEnergy), median(iEnergy), median(cancellation)};
end

function [summary, spectra] = local_schur_audit(matrix, n, directEigenvalues)
eColumns = 1:2*n;
iColumns = 2*n + (1:n);
a = sparse(matrix(eColumns, eColumns));
b = sparse(matrix(eColumns, iColumns));
c = sparse(matrix(iColumns, eColumns));
d = sparse(matrix(iColumns, iColumns));
x = a \ b;
loop = c*x;
h = sparse(d - loop);
y = a \ x;
m = speye(n) + c*y;
hEigenvalues = eig(full(h), 'vector');
generalizedEigenvalues = eig(full(h), full(m), 'vector');

nearestError = zeros(numel(generalizedEigenvalues), 1);
nearestRelativeError = zeros(numel(generalizedEigenvalues), 1);
for index = 1:numel(generalizedEigenvalues)
    [nearestError(index), nearestIndex] = min(abs( ...
        directEigenvalues - generalizedEigenvalues(index)));
    nearestRelativeError(index) = nearestError(index) / ...
        max(abs(directEigenvalues(nearestIndex)), 1e-3);
end

summary = table(condest(a), norm(d, 'fro'), norm(loop, 'fro'), ...
    norm(h, 'fro'), norm(h, 'fro') / (norm(d, 'fro') + norm(loop, 'fro')), ...
    norm(h, 'fro') / norm(d, 'fro'), sprank(h), ...
    svds(h, 1, 'smallest', 'Tolerance', 1e-9, 'MaxIterations', 3000), ...
    sum(abs(hEigenvalues) <= 0.01), sum(abs(hEigenvalues) <= 0.025), ...
    sum(abs(hEigenvalues) <= 0.05), median(nearestError), ...
    prctile(nearestError, 95), median(nearestRelativeError), ...
    prctile(nearestRelativeError, 95), ...
    'VariableNames', {'condestA', 'frobeniusD', 'frobeniusCAinvB', ...
    'frobeniusSchurH', 'schurCancellationRatio', 'schurRelativeToD', ...
    'structuralRankH', 'smallestSingularH', 'countHEigAbsLE1e2', ...
    'countHEigAbsLE2p5e2', 'countHEigAbsLE5e2', ...
    'medianGeneralizedNearestError', 'p95GeneralizedNearestError', ...
    'medianGeneralizedNearestRelativeError', ...
    'p95GeneralizedNearestRelativeError'});
spectra = struct('H', hEigenvalues, 'Generalized', generalizedEigenvalues, ...
    'NearestError', nearestError, 'NearestRelativeError', nearestRelativeError);
end

function definitions = local_ablation_matrices(jRest, j6, jI, n)
j = jRest + j6 + jI;
diagonalBlocks = sparse(size(j, 1), size(j, 2));
for population = 1:3
    indices = (population - 1)*n + (1:n);
    diagonalBlocks(indices, indices) = j(indices, indices);
end
offDiagonalBlocks = j - diagonalBlocks;
eColumns = 1:2*n;
iColumns = 2*n + (1:n);
a = j(eColumns, eColumns);
b = j(eColumns, iColumns);
c = j(iColumns, eColumns);
d = j(iColumns, iColumns);
decoupledEI = blkdiag(a, d);
couplingOnly = [sparse(2*n, 2*n), b; c, sparse(n, n)];
definitions = { ...
    'self_population_blocks_only', diagonalBlocks; ...
    'cross_population_blocks_only', offDiagonalBlocks; ...
    'decoupled_EI_blocks', decoupledEI; ...
    'EI_coupling_only', couplingOnly; ...
    'without_L6', jRest + jI; ...
    'without_inhibition', jRest + j6; ...
    'without_rest', j6 + jI};
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

function row = local_spectrum_row(name, matrix, eigenvalues, n, category)
values = abs(nonzeros(matrix));
if isempty(values); medianEntry = 0; else; medianEntry = median(values); end
row = {name, category, size(matrix, 1), nnz(matrix), ...
    nnz(matrix)/numel(matrix), sprank(matrix), ...
    sum(full(sum(spones(matrix), 2)) == 0), ...
    sum(full(sum(spones(matrix), 1)) == 0), norm(matrix, 'fro'), ...
    full(max(abs(nonzeros(matrix)))), medianEntry, min(abs(eigenvalues)), ...
    max(real(eigenvalues)), min(real(eigenvalues)), max(abs(eigenvalues)), ...
    sum(abs(eigenvalues) <= 1e-10), sum(abs(eigenvalues) <= 1e-4), ...
    sum(abs(eigenvalues) <= 1e-2), sum(abs(eigenvalues) <= 2.5e-2), ...
    sum(abs(eigenvalues) <= 5e-2), sum(abs(real(eigenvalues)) <= 1e-2), ...
    sum(abs(real(eigenvalues)) <= 2.5e-2), ...
    sum(abs(real(eigenvalues)) <= 5e-2)};
if numel(eigenvalues) ~= 3*n
    error('NearZeroCanonical:EigenvalueCount', 'Wrong eigenvalue count for %s.', name);
end
end

function names = local_summary_names()
names = {'matrix', 'category', 'dimension', 'nnz', 'density', ...
    'structuralRank', 'zeroRows', 'zeroColumns', 'frobeniusNorm', ...
    'maximumAbsEntry', 'medianAbsNonzeroEntry', 'minimumAbsEigenvalue', ...
    'maximumRealEigenvalue', 'minimumRealEigenvalue', 'spectralRadius', ...
    'countAbsEigLE1e10', 'countAbsEigLE1e4', 'countAbsEigLE1e2', ...
    'countAbsEigLE2p5e2', 'countAbsEigLE5e2', 'countAbsRealLE1e2', ...
    'countAbsRealLE2p5e2', 'countAbsRealLE5e2'};
end

function eigenvalues = local_block_eigenvalues(matrix, n)
eColumns = 1:2*n;
iColumns = 2*n + (1:n);
tolerance = 1e-13 * max(norm(matrix, 'fro'), 1);
if norm(matrix(:, iColumns), 'fro') <= tolerance
    eigenvalues = [eig(full(matrix(eColumns, eColumns)), 'vector'); zeros(n, 1)];
elseif norm(matrix(:, eColumns), 'fro') <= tolerance
    eigenvalues = [zeros(2*n, 1); ...
        eig(full(matrix(iColumns, iColumns)), 'vector')];
else
    eigenvalues = eig(full(matrix), 'vector');
end
end

function local_write(tableValue, outputRoot, filename)
writetable(tableValue, fullfile(outputRoot, filename), ...
    'FileType', 'text', 'Delimiter', '\t');
end
