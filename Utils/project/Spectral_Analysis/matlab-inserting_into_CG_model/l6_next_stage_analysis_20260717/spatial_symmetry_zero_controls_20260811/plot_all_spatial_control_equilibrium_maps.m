function plot_all_spatial_control_equilibrium_maps(resultRoot)
% Plot S/C/I/E maps for every re-equilibrated spatial manipulation.

if nargin<1 || isempty(resultRoot); resultRoot=getenv('SPATIAL_CONTROL_OUTPUT'); end
outputRoot=fullfile(resultRoot,'equilibrium_maps_all_controls');
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

controls=["baseline";"shared_source_permutation"; ...
    "joint_alpha_0p10";"joint_alpha_0p25";"joint_alpha_0p50"; ...
    "remove_L6_smoothing";"remove_L4_smoothing"; ...
    "remove_L4_and_L6_smoothing";"L4_left_half";"L4_triangle"];
labels=["Baseline";"Shared source permutation"; ...
    "Joint alpha = 0.10";"Joint alpha = 0.25";"Joint alpha = 0.50"; ...
    "Remove L6 smoothing";"Remove L4 smoothing"; ...
    "Remove L4 and L6 smoothing";"L4 one-sided half-Gaussian"; ...
    "L4 directed triangular footprint"];
populations=["S" "C" "I" "E"];
mapSize=[40 40];
states=cell(numel(controls),1);
limits=zeros(numel(populations),2);
limits(:,1)=Inf;
limits(:,2)=-Inf;

for controlIndex=1:numel(controls)
    loaded=load(fullfile(resultRoot,controls(controlIndex)+".mat"),'equilibrium');
    states{controlIndex}=loaded.equilibrium;
    for populationIndex=1:numel(populations)
        values=local_population(loaded.equilibrium,populations(populationIndex));
        limits(populationIndex,1)=min(limits(populationIndex,1),min(values));
        limits(populationIndex,2)=max(limits(populationIndex,2),max(values));
    end
end

for controlIndex=1:numel(controls)
    figureHandle=figure('Color','w','Position',[80 80 1580 430]);
    tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
    for populationIndex=1:numel(populations)
        nexttile
        values=local_population(states{controlIndex},populations(populationIndex));
        imagesc(reshape(values,mapSize));
        axis image off
        colormap(gca,parula);
        clim(limits(populationIndex,:));
        colorbar
        title(populations(populationIndex)+" firing rate",'FontSize',12);
    end
    sgtitle(labels(controlIndex)+" at its re-equilibrated fixed point", ...
        'FontWeight','bold','FontSize',15);
    tag=sprintf('%02d_%s_equilibrium_firing_rate_maps',controlIndex,controls(controlIndex));
    exportgraphics(figureHandle,fullfile(outputRoot,[tag '.pdf']),'ContentType','vector');
    exportgraphics(figureHandle,fullfile(outputRoot,[tag '.png']),'Resolution',180);
    savefig(figureHandle,fullfile(outputRoot,[tag '.fig']));
    close(figureHandle)
end

overview=figure('Color','w','Position',[20 20 1800 2600]);
tiledlayout(numel(controls),numel(populations), ...
    'TileSpacing','compact','Padding','compact');
for controlIndex=1:numel(controls)
    for populationIndex=1:numel(populations)
        nexttile
        values=local_population(states{controlIndex},populations(populationIndex));
        imagesc(reshape(values,mapSize));
        axis image off
        colormap(gca,parula);
        clim(limits(populationIndex,:));
        colorbar
        title(labels(controlIndex)+": "+populations(populationIndex), ...
            'FontSize',8,'Interpreter','none');
    end
end
sgtitle('Firing-rate maps at each manipulation-specific equilibrium', ...
    'FontWeight','bold','FontSize',16);
exportgraphics(overview,fullfile(outputRoot, ...
    '00_all_spatial_manipulations_equilibrium_maps.pdf'),'ContentType','vector');
exportgraphics(overview,fullfile(outputRoot, ...
    '00_all_spatial_manipulations_equilibrium_maps.png'),'Resolution',160);
savefig(overview,fullfile(outputRoot, ...
    '00_all_spatial_manipulations_equilibrium_maps.fig'));
close(overview)

limitsTable=table(populations(:),limits(:,1),limits(:,2), ...
    'VariableNames',{'population','colorMinimum','colorMaximum'});
writetable(limitsTable,fullfile(outputRoot,'shared_color_limits.tsv'), ...
    'FileType','text','Delimiter','\t');
end

function values=local_population(state,population)
switch population
    case "S"
        values=state.S(:);
    case "C"
        values=state.C(:);
    case "I"
        values=state.I(:);
    case "E"
        values=state.S(:)+state.C(:);
end
end
