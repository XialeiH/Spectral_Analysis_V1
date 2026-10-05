function cg = figure1g_prepare_cg_interpolants(R)
% Prepare the condition-independent Paper 3 response interpolants.
l4E = unique(R.l4EMesh);
l4I = unique(R.l4IMesh);
l6 = R.l6Mesh(:);
populations = {'S','C','I'};
cg.F = cell(5,3);
for population = 1:3
    values = R.func.(populations{population});
    for category = 1:5
        cg.F{category,population} = griddedInterpolant( ...
            {l4I,l4E,l6}, squeeze(values(category,:,:,:)), ...
            'linear', 'none');
    end
end
cg.l6Count = numel(l6);
end
