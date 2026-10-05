function run_figure_5D_valley_task
% Compute one point on the Figure 5B minimum-HC regression line.

taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
outputRoot = getenv('FIGURE5D_OUTPUT_ROOT');
setupFile = getenv('FIGURE5D_SETUP_FILE');
if ~isfinite(taskId) || isempty(outputRoot) || ~isfile(setupFile)
    error('Figure5D:Environment', 'Task id, output root, and setup file are required.');
end

regressionSlope = 0.561654445974616;
pathScale = hypot(1, regressionSlope);
anchorBeta6 = [0.073 0.294];
maxArcLength = 0.40 * pathScale;
arcLengthValues = unique([0:0.025:maxArcLength, ...
    anchorBeta6 * pathScale, maxArcLength]);
beta6Values = arcLengthValues / pathScale;
if taskId < 1 || taskId > numel(beta6Values)
    error('Figure5D:TaskId', 'Task id must be in 1:%d.', numel(beta6Values));
end

loaded = load(setupFile, 'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
beta6 = beta6Values(taskId);
betaI = regressionSlope * beta6;
gain6 = 1 + beta6;
gainI = 1 + betaI;
arcLength = beta6 * hypot(1, regressionSlope);

phi = @(state) l6ns_phi(state, 1-gain6, context, [1 1], gainI, 'true');
fixed = real_tuning_fixed_point(phi, baseline, context.RelaxationP);
if ~fixed.Converged
    error('Figure5D:FixedPoint', ...
        'Fixed point failed at beta6 %.6g, betaI %.6g (residual %.3g).', ...
        beta6, betaI, fixed.Residual);
end

n = numel(baseline) / 3;
staticError = HC_norm_diff( ...
    fixed.State(1:n), fixed.State(n+(1:n)), fixed.State(2*n+(1:n)), ...
    baseline(1:n), baseline(n+(1:n)), baseline(2*n+(1:n)), ...
    0.3077, 0.8, 0.2);

J = real_tuning_true_jacobian(fixed.State, context, gain6, gainI);
tauMs = 10.3402405296839;
[maxRealLambda, adaptiveSlowBasis, lambda, vectors] = local_slow_spectrum(J);
spectralEdgePerMs = (maxRealLambda - 1) / tauMs;
if maxRealLambda < 1
    recoveryTimeMs = tauMs / (1 - maxRealLambda);
else
    recoveryTimeMs = NaN;
end

J0 = real_tuning_true_jacobian(baseline, context, 1, 1);
[~, baselineAdaptiveBasis, baselineLambda, baselineVectors] = ...
    local_slow_spectrum(J0);
clusterTolerance = 5e-3;
baselineClusterRank = max(1, sum(real(baselineLambda) >= ...
    real(baselineLambda(1)) - clusterTolerance));
adaptiveClusterRank = max(1, sum(real(lambda) >= ...
    real(lambda(1)) - clusterTolerance));
baselineSlowBasis = orth(baselineVectors(:, 1:baselineClusterRank));
slowBasis = orth(vectors(:, 1:baselineClusterRank));
principalCosines = min(1, max(0, svd(baselineSlowBasis' * slowBasis)));
principalAnglesDeg = acosd(principalCosines);
slowSubspaceAngleDeg = max(principalAnglesDeg);
normalizedProjectorDistance = sqrt(sum(1 - principalCosines.^2) / ...
    baselineClusterRank);
adaptiveSubspaceAngleDeg = rad2deg(subspace( ...
    baselineAdaptiveBasis, adaptiveSlowBasis));
if baselineClusterRank < numel(lambda)
    slowClusterSpectralGap = real(lambda(baselineClusterRank)) - ...
        real(lambda(baselineClusterRank + 1));
else
    slowClusterSpectralGap = NaN;
end

populationRms = [norm(fixed.State(1:n)), ...
    norm(fixed.State(n+(1:n))), norm(fixed.State(2*n+(1:n)))] / sqrt(n);
populationRms = max(populationRms, 1e-8);
weightVector = [ones(n,1)/populationRms(1); ...
    ones(n,1)/populationRms(2); ones(n,1)/populationRms(3)];
cfg = struct('RandomSeed', 7, 'TransientPowerIterations', 18, ...
    'TransientRefinePeak', false);
timeGridMs = unique([0:2.5:20 5:0.25:10 25:5:60 70:10:180]);
transient = l6ns_continuous_transient( ...
    J, timeGridMs, cfg, weightVector, tauMs);

if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end
pointsRoot = fullfile(outputRoot, 'points');
if ~exist(pointsRoot, 'dir'); mkdir(pointsRoot); end
row = table(taskId, beta6, betaI, arcLength, staticError, ...
    maxRealLambda, spectralEdgePerMs, recoveryTimeMs, ...
    transient.PeakGain, transient.PeakTime, slowSubspaceAngleDeg, ...
    normalizedProjectorDistance, adaptiveSubspaceAngleDeg, ...
    baselineClusterRank, adaptiveClusterRank, slowClusterSpectralGap, ...
    fixed.Residual, fixed.Iterations, ...
    'VariableNames', {'taskId','beta6','betaI','arcLength', ...
    'staticResponseError','maxRealLambda','spectralEdgePerMs', ...
    'recoveryTimeMs','peakFiniteTimeGain','peakTimeMs', ...
    'slowSubspaceAngleDeg','normalizedProjectorDistance', ...
    'adaptiveSubspaceAngleDeg','baselineClusterRank', ...
    'adaptiveClusterRank','slowClusterSpectralGap', ...
    'fixedPointResidual','fixedPointIterations'});
outputFile = fullfile(pointsRoot, sprintf('point_%02d.tsv', taskId));
writetable(row, outputFile, 'FileType', 'text', 'Delimiter', '\t');
spectrumTable = table((1:numel(lambda))', real(lambda), imag(lambda), ...
    'VariableNames', {'mode','realLambda','imagLambda'});
writetable(spectrumTable, fullfile(pointsRoot, ...
    sprintf('spectrum_%02d.tsv', taskId)), ...
    'FileType', 'text', 'Delimiter', '\t');
fprintf('Saved %s\n', outputFile);
end

function [maxRealLambda, basis, lambda, vectors] = local_slow_spectrum(J)
options = struct('tol', 1e-9, 'maxit', 1800, 'p', 100, ...
    'isreal', true, 'disp', 0);
try
    [vectors, values] = eigs(J, 20, 'largestreal', options);
catch
    [vectors, values] = eigs(J, 20, 'lr', options);
end
lambda = diag(values);
[~, order] = sort(real(lambda), 'descend');
lambda = lambda(order);
vectors = vectors(:, order);
maxRealLambda = real(lambda(1));
clusterMask = real(lambda) >= maxRealLambda - 5e-3;
basis = orth(vectors(:, clusterMask));
end
