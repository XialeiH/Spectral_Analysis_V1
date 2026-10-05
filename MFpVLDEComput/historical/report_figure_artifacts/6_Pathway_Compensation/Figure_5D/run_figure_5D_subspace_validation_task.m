function run_figure_5D_subspace_validation_task
% Validate slow-subspace rotation without splitting near-degenerate modes.

taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
outputRoot = getenv('FIGURE5D_OUTPUT_ROOT');
setupFile = getenv('FIGURE5D_SETUP_FILE');
if ~isfinite(taskId) || isempty(outputRoot) || ~isfile(setupFile)
    error('Figure5D:Environment', 'Task id, output root, and setup file are required.');
end

beta6Values = [0 0.02 0.04 0.06 0.073 0.08 0.10 0.12 ...
    0.16 0.20 0.24 0.294 0.32 0.36 0.40];
regressionSlope = 0.561654445974616;
beta6 = beta6Values(taskId);
betaI = regressionSlope * beta6;
arcLength = beta6 * hypot(1, regressionSlope);

loaded = load(setupFile, 'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
phi = @(state) l6ns_phi(state, -beta6, context, [1 1], 1+betaI, 'true');
fixed = real_tuning_fixed_point(phi, baseline, context.RelaxationP);
if ~fixed.Converged
    error('Figure5D:FixedPoint', 'Fixed point failed at task %d.', taskId);
end

J = real_tuning_true_jacobian(fixed.State, context, 1+beta6, 1+betaI);
J0 = real_tuning_true_jacobian(baseline, context, 1, 1);
[lambda, vectors] = local_leading_spectrum(J);
[baselineLambda, baselineVectors] = local_leading_spectrum(J0);

ranks = [4 6 8];
angles = zeros(size(ranks));
projectorDistances = zeros(size(ranks));
spectralGaps = zeros(size(ranks));
invariantResiduals = zeros(size(ranks));
principalAnglesRank4Deg = [];
for index = 1:numel(ranks)
    rankUse = ranks(index);
    q0 = orth(baselineVectors(:, 1:rankUse));
    q = orth(vectors(:, 1:rankUse));
    cosines = min(1, max(0, svd(q0' * q)));
    principalAngles = acosd(cosines);
    angles(index) = max(principalAngles);
    if rankUse == 4
        principalAnglesRank4Deg = principalAngles(:)';
    end
    projectorDistances(index) = sqrt(sum(1-cosines.^2) / rankUse);
    spectralGaps(index) = real(lambda(rankUse)) - real(lambda(rankUse+1));
    invariantResiduals(index) = norm(J*q-q*(q'*J*q), 'fro') / ...
        max(norm(J*q, 'fro'), eps);
end

clusterTolerance = 5e-3;
baselineAdaptiveRank = sum(real(baselineLambda) >= ...
    real(baselineLambda(1))-clusterTolerance);
adaptiveRank = sum(real(lambda) >= real(lambda(1))-clusterTolerance);

row = table(taskId, beta6, betaI, arcLength, ...
    angles(1), angles(2), angles(3), ...
    projectorDistances(1), projectorDistances(2), projectorDistances(3), ...
    spectralGaps(1), spectralGaps(2), spectralGaps(3), ...
    invariantResiduals(1), invariantResiduals(2), invariantResiduals(3), ...
    baselineAdaptiveRank, adaptiveRank, fixed.Residual, ...
    'VariableNames', {'taskId','beta6','betaI','arcLength', ...
    'angleRank4Deg','angleRank6Deg','angleRank8Deg', ...
    'projectorRank4','projectorRank6','projectorRank8', ...
    'gapAfterRank4','gapAfterRank6','gapAfterRank8', ...
    'residualRank4','residualRank6','residualRank8', ...
    'baselineAdaptiveRank','adaptiveRank','fixedPointResidual'});
pointsRoot = fullfile(outputRoot, 'points');
if ~exist(pointsRoot, 'dir'); mkdir(pointsRoot); end
writetable(row, fullfile(pointsRoot, sprintf('subspace_%02d.tsv', taskId)), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(pointsRoot, sprintf('subspace_modes_%02d.mat', taskId)), ...
    'taskId', 'beta6', 'betaI', 'arcLength', 'lambda', 'vectors', ...
    'principalAnglesRank4Deg');
end

function [lambda, vectors] = local_leading_spectrum(J)
options = struct('tol', 1e-9, 'maxit', 1800, 'p', 120, ...
    'isreal', true, 'disp', 0);
try
    [vectors, values] = eigs(J, 24, 'largestreal', options);
catch
    [vectors, values] = eigs(J, 24, 'lr', options);
end
lambda = diag(values);
[~, order] = sort(real(lambda), 'descend');
lambda = lambda(order);
vectors = vectors(:, order);
end
