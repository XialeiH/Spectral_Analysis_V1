function aggregate_all_spatial_controls_reequilibrated(outputRoot)
% Aggregate re-equilibrated counts and plot the two footprint-control maps.

if nargin<1 || isempty(outputRoot); outputRoot=getenv('SPATIAL_CONTROL_OUTPUT'); end
controls=["baseline";"shared_source_permutation"; ...
    "joint_alpha_0p10";"joint_alpha_0p25";"joint_alpha_0p50"; ...
    "remove_L6_smoothing";"remove_L4_smoothing"; ...
    "remove_L4_and_L6_smoothing";"L4_left_half";"L4_triangle"];
labels=["baseline";"shared source permutation"; ...
    "joint alpha = 0.10";"joint alpha = 0.25";"joint alpha = 0.50"; ...
    "remove L6 smoothing";"remove L4 smoothing"; ...
    "remove L4 and L6 smoothing";"L4 one-sided half-Gaussian"; ...
    "L4 directed triangular footprint"];
rows=cell(numel(controls),1);
desiredVariables={};
for index=1:numel(controls)
    file=fullfile(outputRoot,controls(index)+".tsv");
    if ~isfile(file); error('Missing result %s.',file); end
    result=readtable(file,'FileType','text','Delimiter','\t');
    if index==1
        desiredVariables=result.Properties.VariableNames;
        desiredVariables{end+1}='solverMethod';
    end
    if ~ismember('solverMethod',result.Properties.VariableNames)
        result.solverMethod=repmat("relaxed_iteration",height(result),1);
    end
    rows{index}=result(:,desiredVariables);
end
summary=vertcat(rows{:});
summary.displayLabel=labels;
writetable(summary,fullfile(outputRoot,'reequilibrated_spatial_control_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

fig=figure('Color','w','Position',[100 100 1450 620]);
axis off
title('Near-zero counts after re-equilibrating every controlled model', ...
    'FontWeight','bold','FontSize',16);
text(0.01,0.90,sprintf('%-42s %8s %8s %8s %12s', ...
    'control','N_{0.01}','N_{0.025}','N_{0.05}','sigma <= 0.05'), ...
    'FontName','Courier','FontWeight','bold','FontSize',12,'Interpreter','none');
for index=1:height(summary)
    line=sprintf('%-42s %8d %8d %8d %12d',labels(index), ...
        summary.countAbsEigLE0p01(index),summary.countAbsEigLE0p025(index), ...
        summary.countAbsEigLE0p05(index),summary.countSingularLE0p05(index));
    text(0.01,0.90-index*0.075,line,'FontName','Courier','FontSize',12, ...
        'Interpreter','none');
end
exportgraphics(fig,fullfile(outputRoot,'12_reequilibrated_near_zero_table.pdf'), ...
    'ContentType','vector');
savefig(fig,fullfile(outputRoot,'12_reequilibrated_near_zero_table.fig'));

baseline=load(fullfile(outputRoot,'baseline.mat'),'equilibrium');
half=load(fullfile(outputRoot,'L4_left_half.mat'),'equilibrium');
triangle=load(fullfile(outputRoot,'L4_triangle.mat'),'equilibrium');
local_plot_maps(outputRoot,baseline.equilibrium,half.equilibrium, ...
    triangle.equilibrium);
end

function local_plot_maps(outputRoot,baseline,half,triangle)
mapSize=[40 40];
populations={'S','C','I','E'};
controlNames={'baseline','one-sided half-Gaussian','directed triangle'};
states={baseline,half,triangle};
fig=figure('Color','w','Position',[40 40 1650 1080]);
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
for row=1:3
    for column=1:4
        nexttile
        values=local_population(states{row},populations{column});
        imagesc(reshape(values,mapSize)); axis image off
        colormap(gca,parula); colorbar
        title(sprintf('%s: %s firing rate',controlNames{row},populations{column}), ...
            'Interpreter','none','FontSize',11);
    end
end
sgtitle('Re-equilibrated firing-rate maps for L4 footprint controls', ...
    'FontWeight','bold');
exportgraphics(fig,fullfile(outputRoot, ...
    '13_L4_footprint_controls_reequilibrated_firing_rate_maps.pdf'), ...
    'ContentType','vector');
exportgraphics(fig,fullfile(outputRoot, ...
    '13_L4_footprint_controls_reequilibrated_firing_rate_maps.png'), ...
    'Resolution',180);
savefig(fig,fullfile(outputRoot, ...
    '13_L4_footprint_controls_reequilibrated_firing_rate_maps.fig'));

fig=figure('Color','w','Position',[40 40 1650 760]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for row=1:2
    state=states{row+1};
    for column=1:4
        nexttile
        difference=local_population(state,populations{column})- ...
            local_population(baseline,populations{column});
        bound=max(abs(difference)); if bound==0; bound=1; end
        imagesc(reshape(difference,mapSize)); axis image off
        colormap(gca,turbo); clim([-bound bound]); colorbar
        title(sprintf('%s - baseline: %s',controlNames{row+1},populations{column}), ...
            'Interpreter','none','FontSize',11);
    end
end
sgtitle('Equilibrium firing-rate changes caused by L4 footprint geometry', ...
    'FontWeight','bold');
exportgraphics(fig,fullfile(outputRoot, ...
    '14_L4_footprint_controls_equilibrium_difference_maps.pdf'), ...
    'ContentType','vector');
exportgraphics(fig,fullfile(outputRoot, ...
    '14_L4_footprint_controls_equilibrium_difference_maps.png'), ...
    'Resolution',180);
savefig(fig,fullfile(outputRoot, ...
    '14_L4_footprint_controls_equilibrium_difference_maps.fig'));
end

function values=local_population(state,name)
if name=="E" || strcmp(name,'E')
    values=state.S(:)+state.C(:);
else
    values=state.(name)(:);
end
end
