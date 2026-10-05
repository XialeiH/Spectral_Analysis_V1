function plot_spatial_shape_reequilibrated_maps(setupFile,resultRoot)
% Compare baseline and re-equilibrated footprint-control firing-rate maps.

loaded = load(setupFile,'setup');
context = loaded.setup.Context;
files = {fullfile(resultRoot,'L4_left_half_reequilibrated.mat'), ...
    fullfile(resultRoot,'L4_triangle_reequilibrated.mat')};
for index=1:numel(files)
    if ~exist(files{index},'file')
        error('SpatialShape:MissingResult','Missing %s.',files{index});
    end
end
left = load(files{1},'equilibrium','summary');
triangle = load(files{2},'equilibrium','summary');

mapSize = context.MapSize;
n = prod(mapSize);
states = {local_unpack(context.FixedPoint(:),n),left.equilibrium,triangle.equilibrium};
rowNames = {'baseline','L4 one-sided half-Gaussian','L4 directed triangular footprint'};
populationNames = {'S','C','I','E=S+C'};
maps = cell(3,4);
for row=1:3
    maps{row,1}=reshape(states{row}.S,mapSize);
    maps{row,2}=reshape(states{row}.C,mapSize);
    maps{row,3}=reshape(states{row}.I,mapSize);
    maps{row,4}=maps{row,1}+maps{row,2};
end

fig = figure('Visible','off','Color','w','Position',[40 40 1750 1200]);
layout = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
for population=1:4
    limits = [min(cellfun(@(x) min(x,[],'all'),maps(:,population))), ...
        max(cellfun(@(x) max(x,[],'all'),maps(:,population)))];
    for row=1:3
        ax=nexttile(layout,(row-1)*4+population);
        imagesc(ax,maps{row,population}); axis(ax,'image'); set(ax,'YDir','normal');
        colormap(ax,parula); colorbar(ax); clim(ax,limits);
        if row==1; title(ax,[populationNames{population} ' firing rate']); end
        if population==1; ylabel(ax,{rowNames{row};'map row'}); end
        xlabel(ax,'map column');
    end
end
title(layout,{sprintf('Re-equilibrated L4 footprint controls: contrast %.0f, orientation %.2f deg', ...
    context.ContrastUse(1),context.OrientationUse(1)), ...
    sprintf('fixed-point residuals: half %.2e; triangle %.2e', ...
    left.summary.fixedPointResidual,triangle.summary.fixedPointResidual)});
base=fullfile(resultRoot,'10_L4_footprint_controls_reequilibrated_firing_rate_maps');
exportgraphics(fig,[base '.png'],'Resolution',220);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
savefig(fig,[base '.fig']); close(fig);

fig = figure('Visible','off','Color','w','Position',[40 40 1750 850]);
layout = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
for population=1:4
    differences={maps{2,population}-maps{1,population}; ...
        maps{3,population}-maps{1,population}};
    maximum=max(cellfun(@(x) max(abs(x),[],'all'),differences));
    for row=1:2
        ax=nexttile(layout,(row-1)*4+population);
        imagesc(ax,differences{row}); axis(ax,'image'); set(ax,'YDir','normal');
        colormap(ax,parula); colorbar(ax); clim(ax,[-maximum maximum]);
        if row==1; title(ax,['Delta ' populationNames{population}]); end
        if population==1; ylabel(ax,{rowNames{row+1};'minus baseline'}); end
        xlabel(ax,'map column');
    end
end
title(layout,'Firing-rate changes caused by re-equilibrating each L4 footprint control');
base=fullfile(resultRoot,'11_L4_footprint_controls_reequilibrated_rate_changes');
exportgraphics(fig,[base '.png'],'Resolution',220);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
savefig(fig,[base '.fig']); close(fig);

save(fullfile(resultRoot,'reequilibrated_firing_rate_maps.mat'), ...
    'maps','rowNames','populationNames','-v7.3');
end

function state=local_unpack(vector,n)
state=struct('S',vector(1:n),'C',vector(n+(1:n)),'I',vector(2*n+(1:n)));
end
