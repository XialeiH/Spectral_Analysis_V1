function result = run_jfull_c_fixed_mode_count_task(taskIndex, existingDataRoot, outputRoot)
% Build a C-population J_full cluster from an exact number of leading modes.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('FIXED_CLUSTER_TASK_INDEX'));
    if ~isfinite(taskIndex)
        taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
    end
end
if nargin < 2 || isempty(existingDataRoot)
    existingDataRoot = getenv('FIXED_CLUSTER_EXISTING_DATA_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('FIXED_CLUSTER_OUTPUT_ROOT');
end

contrasts = [19 42 66 100];
angles = [15 22.5];
if ~isscalar(taskIndex) || ~isfinite(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > numel(contrasts) * numel(angles)
    error('FixedCluster:TaskIndex', 'Task index must be an integer from 1 to 8.');
end
[contrastIndex, angleIndex] = ind2sub([numel(contrasts), numel(angles)], taskIndex);
contrast = contrasts(contrastIndex);
angleDeg = angles(angleIndex);
profile = getenv('FIXED_CLUSTER_PROFILE');
if isempty(profile); profile = 'initial'; end
requestedModeCount = local_mode_count(profile, angleDeg, contrast);
if ~isfolder(existingDataRoot)
    error('FixedCluster:DataRoot', 'Missing existing data root: %s', existingDataRoot);
end
if isempty(outputRoot)
    error('FixedCluster:OutputRoot', 'An output root is required.');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

angleTag = local_angle_tag(angleDeg);
existingFile = fullfile(existingDataRoot, sprintf( ...
    'full_eigenspectrum_maps_angle%s_contrast%d.mat', angleTag, contrast));
existing = load(existingFile, 'conditionData');
sourceFile = existing.conditionData.SourceFiles.DynamicL6;
cacheFile = existing.conditionData.SourceFiles.EigenvalueCache;
if ~isfile(sourceFile) || ~isfile(cacheFile)
    error('FixedCluster:SourceFiles', 'A required source file is missing.');
end

source = load(sourceFile, 'Section4');
cache = load(cacheFile, 'eigenvalueData');
jFull = sparse(source.Section4.A);
cachedEigenvalues = cache.eigenvalueData.J_baseline(:);
dimension = size(jFull, 1);
if ~isequal(size(jFull), [4800 4800]) || numel(cachedEigenvalues) ~= dimension
    error('FixedCluster:Dimension', 'Expected a 4800-dimensional J_full spectrum.');
end

eigsCount = requestedModeCount + 8;
options = struct('tol', 1e-10, 'maxit', 5000, ...
    'p', max(80, 2*eigsCount + 16), 'disp', 0, 'isreal', true);
[vectors, diagonal, flag] = eigs(jFull, eigsCount, 'largestreal', options);
if flag ~= 0
    error('FixedCluster:EigsConvergence', 'eigs returned flag %d.', flag);
end
values = diag(diagonal);
[~, order] = sortrows([real(values), imag(values)], [-1 -2]);
values = values(order);
vectors = vectors(:, order);
cachedTopReal = max(real(cachedEigenvalues));
if abs(real(values(1)) - cachedTopReal) > 1e-8
    error('FixedCluster:TopEigenvalueAudit', ...
        'Recovered top real eigenvalue differs from the cached spectrum.');
end

selectedValues = values(1:requestedModeCount);
selectedVectors = vectors(:, 1:requestedModeCount);
selectedVectors = selectedVectors ./ max(vecnorm(selectedVectors), eps);
[basis, triangularFactor] = qr(selectedVectors, 0);
rankTolerance = max(size(triangularFactor)) * ...
    eps(max(norm(triangularFactor, 2), 1));
clusterRank = sum(abs(diag(triangularFactor)) > rankTolerance);
basis = basis(:, 1:clusterRank);

populationSize = dimension / 3;
mapSide = round(sqrt(populationSize));
cBasis = basis(populationSize + (1:populationSize), :);
cMap = reshape(sqrt(sum(abs(cBasis).^2, 2)), mapSide, mapSide);
residuals = vecnorm(jFull * selectedVectors - ...
    selectedVectors .* reshape(selectedValues, 1, [])) ./ ...
    max((norm(jFull, 'fro') + abs(selectedValues.')) .* ...
    vecnorm(selectedVectors), eps);

result = struct();
result.AngleDeg = angleDeg;
result.Contrast = contrast;
result.Profile = profile;
result.RequestedModeCount = requestedModeCount;
result.ClusterRank = clusterRank;
result.SelectedEigenvalues = selectedValues;
result.CEnvelopeMap = cMap;
result.CutoffRealPartGap = real(values(requestedModeCount)) - ...
    real(values(requestedModeCount + 1));
result.MaximumResidual = max(residuals);
result.SourceFile = sourceFile;
result.EigenvalueCacheFile = cacheFile;

outputFile = fullfile(outputRoot, sprintf( ...
    'Jfull_C_fixed_modes_angle%s_contrast%d.mat', angleTag, contrast));
save(outputFile, 'result', '-v7');
fprintf(['Saved fixed-count J_full C cluster: angle %.2f, contrast %d, ' ...
    'top %d modes, rank %d, cutoff gap %.6g -> %s\n'], ...
    angleDeg, contrast, requestedModeCount, clusterRank, ...
    result.CutoffRealPartGap, outputFile);
end

function tag = local_angle_tag(angleDeg)
tag = strrep(sprintf('%.2f', angleDeg), '.', 'p');
end

function modeCount = local_mode_count(profile, angleDeg, contrast)
switch profile
    case 'initial'
        if angleDeg == 15
            modeCount = 6;
        else
            modeCount = 5;
        end
    case 'revised'
        if angleDeg == 15 && contrast == 66
            modeCount = 10;
        elseif angleDeg == 15 || (angleDeg == 22.5 && ismember(contrast, [19 100]))
            modeCount = 6;
        else
            modeCount = 5;
        end
    otherwise
        error('FixedCluster:Profile', 'Unsupported profile: %s', profile);
end
end
