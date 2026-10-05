function output = figure1g_cg_population_response(recurrentE,recurrentI,l6Axis,weights,interpolants)
% Evaluate the same Paper 3 lookup table without rebuilding its grid.
n = numel(recurrentE);
responses = zeros(n,5);
for category = 1:5
    responses(:,category) = interpolants{category}(recurrentI,recurrentE,l6Axis);
end
output = sum(weights.*responses,2);
end
