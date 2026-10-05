function out = GA_discrete_transient_metrics(A, taskTangents, featureNames, varargin)
% Finite-step transient amplification for the fixed-point update map A.
opts = struct('MaxStep', 20, 'Threshold', 1.05, 'SampleRank', 16, 'Seed', 1);
opts = parse_opts(opts, varargin{:});

n = size(A, 1);
maxStep = opts.MaxStep;
threshold = opts.Threshold;
sampleRank = min(opts.SampleRank, n);
steps = 0:maxStep;

rng(opts.Seed, 'twister');
omega0 = randn(n, sampleRank);
sampleState = omega0;

if nargin < 2 || isempty(taskTangents)
    taskBasis = zeros(n, 0);
    featureState = zeros(n, 0);
else
    taskTangents = full(taskTangents);
    taskBasis = orth(taskTangents);
    featureState = taskTangents;
end
if nargin < 3 || isempty(featureNames)
    featureNames = {};
end

globalGain = nan(numel(steps), 1);
numAmplified = nan(numel(steps), 1);
effectiveDim = nan(numel(steps), 1);
taskGain = nan(numel(steps), 1);
numTaskAmplified = nan(numel(steps), 1);
featureGains = nan(numel(steps), size(featureState, 2));

for si = 1:numel(steps)
    k = steps(si);
    singularVals = randomized_singular_values(A, sampleState, k);
    globalGain(si) = singularVals(1);
    numAmplified(si) = sum(singularVals > threshold);
    ampEnergy = max(singularVals.^2 - 1, 0);
    if any(ampEnergy > 0)
        effectiveDim(si) = sum(ampEnergy)^2 / sum(ampEnergy.^2);
    else
        effectiveDim(si) = 0;
    end

    if ~isempty(taskBasis)
        taskSingularVals = svd(full(taskBasis), 'econ');
        taskGain(si) = taskSingularVals(1);
        numTaskAmplified(si) = sum(taskSingularVals > threshold);
    end
    for fi = 1:size(featureState, 2)
        denom = norm(taskTangents(:, fi));
        if denom > 0
            featureGains(si, fi) = norm(featureState(:, fi)) / denom;
        end
    end

    if si < numel(steps)
        sampleState = A * sampleState;
        if ~isempty(taskBasis)
            taskBasis = A * taskBasis;
        end
        if ~isempty(featureState)
            featureState = A * featureState;
        end
    end
end

[peakGain, peakIdx] = max(globalGain);
[peakTaskGain, peakTaskIdx] = max(taskGain);
[peakFeatureGains, peakFeatureIdx] = max(featureGains, [], 1);

out = struct();
out.Method = 'discrete_fixed_point_power_randomized_range';
out.Propagator = 'A^k';
out.MaxStep = maxStep;
out.TimeGrid = steps;
out.StepGrid = steps;
out.Threshold = threshold;
out.SampleRank = sampleRank;
out.Seed = opts.Seed;
out.GlobalGainByStep = globalGain;
out.PeakGain = peakGain;
out.PeakStep = steps(peakIdx);
out.NumAmplifiedAtPeak = numAmplified(peakIdx);
out.EffectiveDimensionAtPeak = effectiveDim(peakIdx);
out.NumAmplifiedByStep = numAmplified;
out.EffectiveDimensionByStep = effectiveDim;
out.TaskGainByStep = taskGain;
out.PeakTaskGain = peakTaskGain;
out.PeakTaskStep = steps(peakTaskIdx);
out.NumTaskAmplifiedAtTaskPeak = numTaskAmplified(peakTaskIdx);
out.NumTaskAmplifiedAtGlobalPeak = numTaskAmplified(peakIdx);
out.FeatureNames = featureNames;
out.FeatureGainByStep = featureGains;
out.PeakFeatureGains = peakFeatureGains;
out.PeakFeatureSteps = steps(peakFeatureIdx);
end

function s = randomized_singular_values(A, rangeState, k)
[Q, ~] = qr(full(rangeState), 0);
adjointState = Q;
for ii = 1:k
    adjointState = A' * adjointState;
end
s = svd(full(adjointState), 'econ');
if isempty(s)
    s = NaN;
end
end

function opts = parse_opts(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    value = varargin{k+1};
    if ~isfield(opts, name)
        error('Unknown option %s.', name);
    end
    opts.(name) = value;
end
end
