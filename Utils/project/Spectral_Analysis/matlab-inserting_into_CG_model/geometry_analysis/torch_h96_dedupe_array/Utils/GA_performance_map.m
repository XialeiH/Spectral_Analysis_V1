function perf = GA_performance_map(section4, varargin)
% Section-9 scalar performance summaries assembled from geometry outputs.
opts = struct('TaskMetric', [], 'HCNorm', [], 'Extra', struct());
opts = parse_opts(opts, varargin{:});

perf.MaxRealLambda = section4.MaxRealLambda;
perf.StabilityMargin = section4.StabilityMargin;
perf.SpectralRadius = section4.SpectralRadius;
perf.NonNormality = section4.NonNormality;

if ~isempty(opts.TaskMetric)
    metricSummary = GA_task_metric_summary(opts.TaskMetric);
    perf.TaskMetricTrace = metricSummary.Trace;
    perf.TaskMetricConditionNumber = metricSummary.ConditionNumber;
    perf.TaskMetricMinEigenValue = metricSummary.MinEigenValue;
else
    perf.TaskMetricTrace = [];
    perf.TaskMetricConditionNumber = [];
    perf.TaskMetricMinEigenValue = [];
end

perf.HCNorm = opts.HCNorm;
perf.Extra = opts.Extra;
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
