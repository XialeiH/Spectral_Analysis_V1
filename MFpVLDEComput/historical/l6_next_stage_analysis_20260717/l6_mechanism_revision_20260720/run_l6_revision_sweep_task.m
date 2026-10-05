function summary = run_l6_revision_sweep_task(taskIndex, sourceRoot, outputRoot)
% Reconstruct the affine L6-weight Jacobian and sweep its leading eigenvalue.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('L6R_SOURCE_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('L6R_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > numel(angles) * numel(contrasts)
    error('L6Revision:SweepTask', 'Task index must be an integer from 1 to 16.');
end
if isempty(sourceRoot) || ~isfolder(sourceRoot)
    error('L6Revision:SourceRoot', 'L6R_SOURCE_ROOT must identify the mechanism run root.');
end
if isempty(outputRoot)
    error('L6Revision:OutputRoot', 'L6R_OUTPUT_ROOT is required.');
end

angleIndex = floor((taskIndex - 1) / numel(contrasts)) + 1;
contrastIndex = mod(taskIndex - 1, numel(contrasts)) + 1;
angle = angles(angleIndex);
contrast = contrasts(contrastIndex);

[jacobian0, jacobian1] = local_endpoint_jacobians(sourceRoot, angle, contrast);
jacobianDelta = jacobian1 - jacobian0;

w = (-0.30:0.025:1.50)';
boundaryFile = getenv('L6R_BOUNDARY_TSV');
if ~isempty(boundaryFile) && isfile(boundaryFile)
    boundary = readtable(boundaryFile, 'FileType', 'text', 'Delimiter', '\t');
    row = boundary.angle == angle & boundary.contrast == contrast;
    if any(row)
        w = unique([w; boundary.lowW(row); boundary.w_boundary(row); boundary.highW(row)]);
    end
end
w = sort(w);

maxRealLambda = nan(size(w));
stabilityMargin = nan(size(w));
eigsFlag = nan(size(w));
runtimeSeconds = nan(size(w));
opts = struct('tol', 1e-9, 'maxit', 1800, 'p', 64, 'disp', 0);
warmStart = [];

for index = 1:numel(w)
    jacobian = jacobian0 + w(index) * jacobianDelta;
    if ~isempty(warmStart)
        opts.v0 = warmStart;
    end
    timer = tic;
    [vectors, values, flag] = eigs(jacobian, 6, 'largestreal', opts);
    runtimeSeconds(index) = toc(timer);
    eigenvalues = diag(values);
    [maxRealLambda(index), leadingIndex] = max(real(eigenvalues));
    stabilityMargin(index) = 1 - maxRealLambda(index);
    eigsFlag(index) = flag;
    warmStart = vectors(:, leadingIndex) / norm(vectors(:, leadingIndex));
    fprintf('angle %.2f contrast %d w %+.9f: max Re(lambda) %.12g, flag %d.\n', ...
        angle, contrast, w(index), maxRealLambda(index), flag);
end

if any(~isfinite(maxRealLambda)) || any(eigsFlag ~= 0)
    error('L6Revision:EigsFailure', 'At least one leading-eigenvalue solve failed.');
end

summary = table(repmat(angle, size(w)), repmat(contrast, size(w)), w, 1-w, ...
    maxRealLambda, stabilityMargin, eigsFlag, runtimeSeconds, ...
    'VariableNames', {'angle','contrast','w','gammaL6','maxRealLambda', ...
    'stabilityMargin','eigsFlag','runtimeSeconds'});

dataRoot = fullfile(outputRoot, 'sweep_data');
if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end
stem = sprintf('l6_weight_sweep_angle%.2f_contrast%d', angle, contrast);
writetable(summary, fullfile(dataRoot, [stem '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
metadata = struct('Definition', 'J(w)=J(0)+w*(J(1)-J(0))', ...
    'W0Meaning', 'original dynamic L6', 'W1Meaning', 'frozen L6', ...
    'SourceRoot', sourceRoot, 'TaskIndex', taskIndex);
save(fullfile(dataRoot, [stem '.mat']), 'summary', 'metadata');
end

function [jacobian0, jacobian1] = local_endpoint_jacobians(sourceRoot, angle, contrast)
file0 = local_endpoint_file(sourceRoot, angle, contrast, 0);
file1 = local_endpoint_file(sourceRoot, angle, contrast, 1);
loaded0 = load(file0, 'Section4');
loaded1 = load(file1, 'Section4');
jacobian0 = loaded0.Section4.A;
jacobian1 = loaded1.Section4.A;
if ~isequal(size(jacobian0), size(jacobian1)) || size(jacobian0,1) ~= size(jacobian0,2)
    error('L6Revision:JacobianSize', 'Endpoint Jacobians have incompatible dimensions.');
end
end

function file = local_endpoint_file(sourceRoot, angle, contrast, w)
tag = sprintf('L6eqW%.2f', w);
tag = strrep(strrep(tag, '.', 'p'), '-', 'm');
file = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_%s_contr%d_angle_%.2f.mat', ...
    tag, contrast, angle));
if ~isfile(file)
    error('L6Revision:MissingEndpoint', 'Missing endpoint file %s.', file);
end
end
