function outputs = run_near_zero_condition_task(taskIndex, sourceRoot, outputRoot)
% Audit near-zero spectra for the two-pathway Jacobian decomposition.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('NZ_SOURCE_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('NZ_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
[contrastIndex, angleIndex] = ind2sub([numel(contrasts), numel(angles)], taskIndex);
angleValue = angles(angleIndex);
contrast = contrasts(contrastIndex);

if taskIndex < 1 || taskIndex > numel(angles) * numel(contrasts) || ...
        taskIndex ~= round(taskIndex)
    error('NearZero:TaskIndex', 'Task index must be an integer from 1 to 16.');
end
if ~isfolder(sourceRoot)
    error('NearZero:SourceRoot', 'Source root does not exist: %s', sourceRoot);
end
if isempty(outputRoot)
    error('NearZero:OutputRoot', 'Output root is required.');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

file0 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW0p00_contr%d_angle_%.2f.mat', ...
    contrast, angleValue));
file1 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW1p00_contr%d_angle_%.2f.mat', ...
    contrast, angleValue));
local_require_file(file0);
local_require_file(file1);

data0 = load(file0, 'Section4');
data1 = load(file1, 'Section4');
jBaseline = sparse(data0.Section4.A);
jNoL6 = sparse(data1.Section4.A);
dimension = size(jBaseline, 1);
if dimension ~= 4800 || size(jBaseline, 2) ~= dimension
    error('NearZero:Dimension', 'Expected a 4800-by-4800 Jacobian.');
end
populationSize = dimension / 3;
iColumns = 2 * populationSize + (1:populationSize);

j6 = jBaseline - jNoL6;
jI = sparse(dimension, dimension);
jI(:, iColumns) = jBaseline(:, iColumns);
jRest = jBaseline - j6 - jI;

reconstructionError = norm(jRest + j6 + jI - jBaseline, 'fro') / ...
    max(norm(jBaseline, 'fro'), eps);
l6ISourceFraction = norm(j6(:, iColumns), 'fro') / max(norm(j6, 'fro'), eps);
iNonISourceFraction = norm(jI(:, 1:2*populationSize), 'fro') / ...
    max(norm(jI, 'fro'), eps);
if reconstructionError > 1e-12 || l6ISourceFraction > 1e-12 || ...
        iNonISourceFraction > 1e-12
    error('NearZero:Decomposition', ...
        'Two-pathway decomposition failed its source-column audit.');
end

matrices = { ...
    'J_baseline', jBaseline, local_saved_eigenvalues(data0.Section4); ...
    'J_rest', jRest, []; ...
    'J_6', j6, []; ...
    'J_I', jI, []; ...
    'J_rest_plus_J6', jRest + j6, []; ...
    'J_rest_plus_JI', jNoL6, local_saved_eigenvalues(data1.Section4); ...
    'J6_plus_JI', j6 + jI, []};

matrixRows = cell(size(matrices, 1), 31);
blockRows = cell(0, 9);
eigenvalueData = struct();
for matrixIndex = 1:size(matrices, 1)
    name = matrices{matrixIndex, 1};
    matrix = sparse(matrices{matrixIndex, 2});
    eigenvalues = matrices{matrixIndex, 3};
    if isempty(eigenvalues)
        eigenvalues = local_block_eigenvalues(matrix, populationSize);
    end
    eigenvalues = eigenvalues(:);
    eigenvalueData.(name) = eigenvalues;
    matrixRows(matrixIndex, :) = local_matrix_row(name, matrix, eigenvalues, ...
        angleValue, contrast, populationSize);
    blockRows = [blockRows; local_block_rows(name, matrix, angleValue, ...
        contrast, populationSize)]; %#ok<AGROW>
end

matrixSummary = cell2table(matrixRows, 'VariableNames', { ...
    'matrix', 'angleDeg', 'contrast', 'dimension', 'nnz', 'density', ...
    'structuralRank', 'zeroRows', 'zeroColumns', 'diagonalNnz', ...
    'frobeniusNorm', 'oneNorm', 'infinityNorm', 'maximumAbsEntry', ...
    'medianAbsNonzeroEntry', 'minimumAbsEigenvalue', 'maximumRealEigenvalue', ...
    'minimumRealEigenvalue', 'spectralRadius', 'countAbsEigLE1e10', ...
    'countAbsEigLE1e6', 'countAbsEigLE1e4', 'countAbsEigLE1e3', ...
    'countAbsEigLE1e2', 'countAbsEigLE2p5e2', 'countAbsEigLE5e2', ...
    'countAbsRealLE1e2', 'countAbsRealLE2p5e2', 'countAbsRealLE5e2', ...
    'fractionAbsEigLE5e2', 'fractionAbsRealLE5e2'});
blockSummary = cell2table(blockRows, 'VariableNames', { ...
    'matrix', 'angleDeg', 'contrast', 'targetPopulation', 'sourcePopulation', ...
    'nnz', 'density', 'frobeniusNorm', 'maximumAbsEntry'});

gamma6Values = [-3 -1 0 0.5 1 2 4];
proportionalityRows = cell(numel(gamma6Values) + 3, 9);
rowIndex = 0;
for gamma6 = gamma6Values
    rowIndex = rowIndex + 1;
    matrix = jRest + gamma6 * j6;
    proportionalityRows(rowIndex, :) = local_proportionality_row( ...
        sprintf('J_rest_plus_%+.3g_J6', gamma6), matrix, angleValue, ...
        contrast, gamma6, populationSize);
end
extraMatrices = {'J_rest', jRest; 'J_6', j6; 'J_baseline', jBaseline};
for entryIndex = 1:size(extraMatrices, 1)
    rowIndex = rowIndex + 1;
    proportionalityRows(rowIndex, :) = local_proportionality_row( ...
        extraMatrices{entryIndex, 1}, extraMatrices{entryIndex, 2}, ...
        angleValue, contrast, NaN, populationSize);
end
proportionalitySummary = cell2table(proportionalityRows, 'VariableNames', { ...
    'matrix', 'angleDeg', 'contrast', 'gamma6', 'bestScalarQ', ...
    'globalRelativeResidual', 'medianColumnScalarQ', ...
    'minimumColumnScalarQ', 'maximumColumnScalarQ'});

audit = table(angleValue, contrast, reconstructionError, l6ISourceFraction, ...
    iNonISourceFraction, 'VariableNames', {'angleDeg', 'contrast', ...
    'reconstructionRelativeError', 'J6ISourceFraction', ...
    'JINonISourceFraction'});

tag = sprintf('angle%s_contrast%d', local_number_tag(angleValue), contrast);
writetable(matrixSummary, fullfile(outputRoot, ...
    ['near_zero_matrix_summary_', tag, '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
writetable(blockSummary, fullfile(outputRoot, ...
    ['near_zero_block_summary_', tag, '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
writetable(proportionalitySummary, fullfile(outputRoot, ...
    ['near_zero_proportionality_', tag, '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
writetable(audit, fullfile(outputRoot, ['near_zero_audit_', tag, '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(outputRoot, ['near_zero_eigenvalues_', tag, '.mat']), ...
    'eigenvalueData', 'angleValue', 'contrast', '-v7.3');

outputs = struct('MatrixSummary', matrixSummary, 'BlockSummary', blockSummary, ...
    'ProportionalitySummary', proportionalitySummary, 'Audit', audit);
fprintf('Near-zero condition audit complete: angle %.2f, contrast %d.\n', ...
    angleValue, contrast);
end

function row = local_matrix_row(name, matrix, eigenvalues, angleValue, ...
        contrast, populationSize)
dimension = size(matrix, 1);
values = abs(nonzeros(matrix));
if isempty(values); medianEntry = 0; else; medianEntry = median(values); end
row = {name, angleValue, contrast, dimension, nnz(matrix), ...
    nnz(matrix) / numel(matrix), sprank(matrix), ...
    sum(full(sum(spones(matrix), 2)) == 0), ...
    sum(full(sum(spones(matrix), 1)) == 0), nnz(diag(matrix)), ...
    norm(matrix, 'fro'), norm(matrix, 1), norm(matrix, inf), ...
    full(max(abs(nonzeros(matrix)))), medianEntry, min(abs(eigenvalues)), ...
    max(real(eigenvalues)), min(real(eigenvalues)), max(abs(eigenvalues)), ...
    sum(abs(eigenvalues) <= 1e-10), sum(abs(eigenvalues) <= 1e-6), ...
    sum(abs(eigenvalues) <= 1e-4), sum(abs(eigenvalues) <= 1e-3), ...
    sum(abs(eigenvalues) <= 1e-2), sum(abs(eigenvalues) <= 2.5e-2), ...
    sum(abs(eigenvalues) <= 5e-2), sum(abs(real(eigenvalues)) <= 1e-2), ...
    sum(abs(real(eigenvalues)) <= 2.5e-2), ...
    sum(abs(real(eigenvalues)) <= 5e-2), ...
    sum(abs(eigenvalues) <= 5e-2) / (3 * populationSize), ...
    sum(abs(real(eigenvalues)) <= 5e-2) / (3 * populationSize)};
end

function rows = local_block_rows(name, matrix, angleValue, contrast, n)
labels = {'S', 'C', 'I'};
rows = cell(9, 9);
rowIndex = 0;
for target = 1:3
    targetRows = (target - 1) * n + (1:n);
    for source = 1:3
        rowIndex = rowIndex + 1;
        sourceColumns = (source - 1) * n + (1:n);
        block = matrix(targetRows, sourceColumns);
        entries = abs(nonzeros(block));
        if isempty(entries); maximumEntry = 0; else; maximumEntry = max(entries); end
        rows(rowIndex, :) = {name, angleValue, contrast, labels{target}, ...
            labels{source}, nnz(block), nnz(block) / numel(block), ...
            norm(block, 'fro'), maximumEntry};
    end
end
end

function row = local_proportionality_row(name, matrix, angleValue, ...
        contrast, gamma6, n)
sColumns = matrix(:, 1:n);
cColumns = matrix(:, n + (1:n));
denominator = full(sum(sum(sColumns .* sColumns)));
if denominator <= eps
    q = NaN;
    residual = NaN;
    columnQ = nan(n, 1);
else
    q = full(sum(sum(sColumns .* cColumns))) / denominator;
    residual = norm(cColumns - q * sColumns, 'fro') / ...
        max(norm(cColumns, 'fro'), eps);
    columnQ = nan(n, 1);
    for column = 1:n
        source = sColumns(:, column);
        target = cColumns(:, column);
        sourceNormSquared = full(source' * source);
        if sourceNormSquared > eps
            columnQ(column) = full(source' * target) / sourceNormSquared;
        end
    end
end
valid = columnQ(isfinite(columnQ));
if isempty(valid)
    medianQ = NaN; minimumQ = NaN; maximumQ = NaN;
else
    medianQ = median(valid); minimumQ = min(valid); maximumQ = max(valid);
end
row = {name, angleValue, contrast, gamma6, q, residual, medianQ, ...
    minimumQ, maximumQ};
end

function eigenvalues = local_block_eigenvalues(matrix, n)
dimension = size(matrix, 1);
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
if numel(eigenvalues) ~= dimension
    error('NearZero:EigenvalueCount', 'Expected %d eigenvalues.', dimension);
end
end

function eigenvalues = local_saved_eigenvalues(section4)
eigenvalues = [];
if isfield(section4, 'EigenValues') && numel(section4.EigenValues) == 4800
    eigenvalues = section4.EigenValues(:);
end
end

function tag = local_number_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
tag = strrep(tag, '-', 'm');
end

function local_require_file(path)
if ~isfile(path)
    error('NearZero:MissingFile', 'Required file does not exist: %s', path);
end
end
