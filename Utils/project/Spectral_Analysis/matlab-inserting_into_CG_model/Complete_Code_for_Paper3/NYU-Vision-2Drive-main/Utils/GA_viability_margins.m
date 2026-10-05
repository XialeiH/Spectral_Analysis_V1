function margins = GA_viability_margins(perf, req)
% Convert requirement bounds into signed margins. Positive means viable.
names = fieldnames(req);
for k = 1:numel(names)
    name = names{k};
    bound = req.(name);
    value = perf.(name);
    if isfield(bound, 'Min')
        margins.(name) = value - bound.Min;
    elseif isfield(bound, 'Max')
        margins.(name) = bound.Max - value;
    else
        error('Requirement %s must have Min or Max.', name);
    end
end
margVals = zeros(1, numel(names));
for k = 1:numel(names)
    margVals(k) = margins.(names{k});
end
margins.AllMargins = margVals;
margins.MinMargin = min(margVals);
margins.IsViable = all(margVals >= 0);
end
