function plot_spatial_shape_control_firing_rate_maps(setupFile,outputRoot)
% Plot the frozen canonical rates used by the two L4 footprint controls.

loaded = load(setupFile,'setup');
context = loaded.setup.Context;
state = context.FixedPoint(:);
mapSize = context.MapSize;
n = prod(mapSize);
if numel(state) ~= 3*n
    error('SpatialShape:FixedPoint','Expected a three-population fixed point.');
end

S = reshape(state(1:n),mapSize);
C = reshape(state(n+(1:n)),mapSize);
I = reshape(state(2*n+(1:n)),mapSize);
E = S+C;
maps = {S,C,I,E};
populationNames = {'S','C','I','E=S+C'};
controlNames = {'L4 one-sided half-Gaussian','L4 directed triangular footprint'};

fig = figure('Visible','off','Color','w','Position',[50 50 1700 850]);
layout = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
for row = 1:2
    for population = 1:4
        ax = nexttile(layout);
        imagesc(ax,maps{population});
        axis(ax,'image');
        set(ax,'YDir','normal');
        colormap(ax,parula);
        colorbar(ax);
        clim(ax,[min(maps{population},[],'all'),max(maps{population},[],'all')]);
        if row==1
            title(ax,sprintf('%s firing rate',populationNames{population}));
        end
        if population==1
            ylabel(ax,{controlNames{row};'map row'});
        end
        xlabel(ax,'map column');
    end
end
title(layout,{sprintf('Canonical firing-rate maps: contrast %.0f, orientation %.2f deg', ...
    context.ContrastUse(1),context.OrientationUse(1)), ...
    'Rows are identical by construction: footprint controls modify only the fixed-point Jacobian derivatives'});

if ~exist(outputRoot,'dir'); mkdir(outputRoot); end
base = fullfile(outputRoot,'09_L4_footprint_controls_frozen_firing_rate_maps');
exportgraphics(fig,[base '.png'],'Resolution',220);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
savefig(fig,[base '.fig']);
close(fig);

metadata = struct('Contrast',context.ContrastUse(1), ...
    'Orientation',context.OrientationUse(1), ...
    'Controls',{controlNames},'RowsIdenticalByConstruction',true, ...
    'Intervention','fixed-operating-point Jacobian derivative only');
save([base '.mat'],'S','C','I','E','metadata');
end
