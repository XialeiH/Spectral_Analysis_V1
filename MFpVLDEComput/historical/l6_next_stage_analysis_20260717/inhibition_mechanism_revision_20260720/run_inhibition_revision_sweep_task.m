function [summary, boundary] = run_inhibition_revision_sweep_task(taskIndex, sourceRoot, outputRoot)
% Sweep the affine inhibition-freeze Jacobian and locate its stability boundary.

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
contrasts = [19 42 66 100];
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > numel(angles) * numel(contrasts)
    error('InhibitionRevision:SweepTask', ...
        'Task index must be an integer from 1 to 16.');
end
if isempty(sourceRoot) || ~isfolder(sourceRoot)
    error('InhibitionRevision:SourceRoot', ...
        'INHR_SOURCE_ROOT must identify the celltype-control run root.');
end
if isempty(outputRoot)
    error('InhibitionRevision:OutputRoot', 'INHR_OUTPUT_ROOT is required.');
end

angleIndex = floor((taskIndex - 1) / numel(contrasts)) + 1;
contrastIndex = mod(taskIndex - 1, numel(contrasts)) + 1;
angle = angles(angleIndex);
contrast = contrasts(contrastIndex);

[jacobian0, jacobian1] = local_endpoint_jacobians(sourceRoot, angle, contrast);
jacobianDelta = jacobian1 - jacobian0;

wGrid = (-0.30:0.025:1.50)';
maxRealLambda = nan(size(wGrid));
eigsFlag = nan(size(wGrid));
runtimeSeconds = nan(size(wGrid));
opts = struct('tol', 1e-9, 'maxit', 1800, 'p', 64, 'disp', 0);
warmStart = [];
for index = 1:numel(wGrid)
    [maxRealLambda(index), eigsFlag(index), warmStart, runtimeSeconds(index)] = ...
        local_leading_eigenvalue(jacobian0 + wGrid(index) * jacobianDelta, ...
        opts, warmStart);
    fprintf(['angle %.2f contrast %d w_I %+.9f: max Re(lambda) ' ...
        '%.12g, flag %d.\n'], angle, contrast, wGrid(index), ...
        maxRealLambda(index), eigsFlag(index));
end
if any(~isfinite(maxRealLambda)) || any(eigsFlag ~= 0)
    error('InhibitionRevision:EigsFailure', ...
        'At least one broad-sweep leading-eigenvalue solve failed.');
end

margin = 1 - maxRealLambda;
crossingIndex = find(margin(1:end-1) >= 0 & margin(2:end) <= 0, 1, 'first');
if isempty(crossingIndex)
    error('InhibitionRevision:MissingCrossing', ...
        ['No stable-to-unstable crossing was found over w_I in [-0.30, 1.50] ' ...
        'for angle %.2f, contrast %d.'], angle, contrast);
end

lowW = wGrid(crossingIndex);
lowMaxReal = maxRealLambda(crossingIndex);
highW = wGrid(crossingIndex + 1);
highMaxReal = maxRealLambda(crossingIndex + 1);
boundaryIterations = 0;
boundaryWarmStart = warmStart;
while highW - lowW > 1e-6 && boundaryIterations < 24
    boundaryIterations = boundaryIterations + 1;
    midpoint = 0.5 * (lowW + highW);
    [midMaxReal, midFlag, boundaryWarmStart] = local_leading_eigenvalue( ...
        jacobian0 + midpoint * jacobianDelta, opts, boundaryWarmStart);
    if midFlag ~= 0 || ~isfinite(midMaxReal)
        error('InhibitionRevision:BoundaryEigsFailure', ...
            'Boundary eigensolve failed at w_I %.12g.', midpoint);
    end
    if midMaxReal < 1
        lowW = midpoint;
        lowMaxReal = midMaxReal;
    else
        highW = midpoint;
        highMaxReal = midMaxReal;
    end
end
wBoundary = 0.5 * (lowW + highW);
[boundaryMaxReal, boundaryFlag, ~, boundaryRuntime] = local_leading_eigenvalue( ...
    jacobian0 + wBoundary * jacobianDelta, opts, boundaryWarmStart);
if boundaryFlag ~= 0 || ~isfinite(boundaryMaxReal)
    error('InhibitionRevision:BoundaryEigsFailure', ...
        'Final boundary eigensolve failed at w_I %.12g.', wBoundary);
end
if ~(lowMaxReal < 1 && highMaxReal >= 1)
    error('InhibitionRevision:InvalidBracket', ...
        'The final stability bracket does not straddle max Re(lambda)=1.');
end

extraW = [lowW; wBoundary; highW];
extraMaxReal = [lowMaxReal; boundaryMaxReal; highMaxReal];
extraFlag = zeros(3,1);
extraRuntime = [nan; boundaryRuntime; nan];
w = [wGrid; extraW];
maxRealLambda = [maxRealLambda; extraMaxReal];
eigsFlag = [eigsFlag; extraFlag];
runtimeSeconds = [runtimeSeconds; extraRuntime];
[w, order] = sort(w);
maxRealLambda = maxRealLambda(order);
eigsFlag = eigsFlag(order);
runtimeSeconds = runtimeSeconds(order);
[~, uniqueIndex] = unique(w, 'stable');
w = w(uniqueIndex);
maxRealLambda = maxRealLambda(uniqueIndex);
eigsFlag = eigsFlag(uniqueIndex);
runtimeSeconds = runtimeSeconds(uniqueIndex);

summary = table(repmat(angle, size(w)), repmat(contrast, size(w)), w, 1-w, ...
    maxRealLambda, 1-maxRealLambda, eigsFlag, runtimeSeconds, ...
    'VariableNames', {'angle','contrast','wI','gammaI','maxRealLambda', ...
    'stabilityMargin','eigsFlag','runtimeSeconds'});
boundary = table(angle, contrast, wBoundary, 1-wBoundary, boundaryMaxReal, ...
    lowW, lowMaxReal, highW, highMaxReal, boundaryIterations, highW-lowW, ...
    'VariableNames', {'angle','contrast','wBoundary','gammaBoundary', ...
    'maxRealAtBoundary','lowW','lowMaxReal','highW','highMaxReal', ...
    'iterations','bracketWidth'});

sweepRoot = fullfile(outputRoot, 'sweep_data');
boundaryRoot = fullfile(outputRoot, 'boundary_data');
if ~exist(sweepRoot, 'dir'); mkdir(sweepRoot); end
if ~exist(boundaryRoot, 'dir'); mkdir(boundaryRoot); end
sweepStem = sprintf('inhibition_weight_sweep_angle%.2f_contrast%d', ...
    angle, contrast);
boundaryStem = sprintf('inhibition_boundary_angle%.2f_contrast%d', ...
    angle, contrast);
writetable(summary, fullfile(sweepRoot, [sweepStem '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
writetable(boundary, fullfile(boundaryRoot, [boundaryStem '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
metadata = struct('Definition', 'J_I(w_I)=J_I(0)+w_I*(J_I(1)-J_I(0))', ...
    'W0Meaning', 'original dynamic inhibition', ...
    'W1Meaning', 'frozen inhibitory derivative with dynamic L6', ...
    'SourceRoot', sourceRoot, 'TaskIndex', taskIndex);
save(fullfile(sweepRoot, [sweepStem '.mat']), 'summary', 'boundary', 'metadata');
fprintf(['Boundary angle %.2f contrast %d: stable w_I %.12g ' ...
    '(%.12g), unstable w_I %.12g (%.12g).\n'], angle, contrast, ...
    lowW, lowMaxReal, highW, highMaxReal);
end

function [maxReal, flag, leadingVector, runtimeSeconds] = ...
        local_leading_eigenvalue(jacobian, opts, warmStart)
if nargin >= 3 && ~isempty(warmStart)
    opts.v0 = warmStart;
elseif isfield(opts, 'v0')
    opts = rmfield(opts, 'v0');
end
timer = tic;
[vectors, values, flag] = eigs(jacobian, 6, 'largestreal', opts);
runtimeSeconds = toc(timer);
eigenvalues = diag(values);
[maxReal, leadingIndex] = max(real(eigenvalues));
leadingVector = vectors(:, leadingIndex);
leadingVector = leadingVector / max(norm(leadingVector), eps);
end

function [jacobian0, jacobian1] = local_endpoint_jacobians(sourceRoot, angle, contrast)
file0 = local_endpoint_file(sourceRoot, angle, contrast, 0);
file1 = local_endpoint_file(sourceRoot, angle, contrast, 1);
loaded0 = load(file0, 'Section4');
loaded1 = load(file1, 'Section4');
jacobian0 = loaded0.Section4.A;
jacobian1 = loaded1.Section4.A;
if ~isequal(size(jacobian0), size(jacobian1)) || ...
        size(jacobian0,1) ~= size(jacobian0,2)
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
