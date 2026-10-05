function summary = GA_task_metric_summary(g)
% Basic section-6 summaries for a task-space metric.
g = 0.5 * (g + g');
ev = eig(full(g));
ev = sort(real(ev), 'descend');

summary.Metric = g;
summary.EigenValues = ev;
summary.Trace = trace(g);
summary.Det = det(full(g));
summary.MinEigenValue = min(ev);
summary.MaxEigenValue = max(ev);
summary.ConditionNumber = summary.MaxEigenValue / max(summary.MinEigenValue, eps);

d = sqrt(max(diag(g), eps));
summary.Correlation = g ./ (d * d');
summary.OffDiagonalNorm = norm(summary.Correlation - diag(diag(summary.Correlation)), 'fro');
end
