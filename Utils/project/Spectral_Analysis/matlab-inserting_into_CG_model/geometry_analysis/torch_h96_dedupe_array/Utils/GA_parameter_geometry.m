function out = GA_parameter_geometry(Dur, Dpr, Sigma, varargin)
% Section-8 parameter geometry blocks H and C with the task block g.
out = GA_reliability_geometry(Dur, Sigma, 'Dpr', Dpr, varargin{:});
out.TaskMetricSummary = GA_task_metric_summary(out.g);
out.ParameterMetricSummary = GA_task_metric_summary(out.H);
out.ConfoundNorm = norm(out.C, 'fro');
end
