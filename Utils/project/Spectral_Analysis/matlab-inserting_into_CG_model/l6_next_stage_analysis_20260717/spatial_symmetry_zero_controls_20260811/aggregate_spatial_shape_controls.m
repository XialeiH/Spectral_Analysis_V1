function summary = aggregate_spatial_shape_controls(outputRoot)
% Aggregate and plot circular/left-half/triangle spatial controls.

if nargin < 1 || isempty(outputRoot)
    outputRoot = getenv('SPATIAL_CONTROL_OUTPUT');
end
files = dir(fullfile(outputRoot,'shape_*.tsv'));
if numel(files)~=7
    error('SpatialShape:Missing','Expected 7 summaries, found %d.',numel(files));
end
parts = cell(numel(files),1);
for index=1:numel(files)
    parts{index}=readtable(fullfile(files(index).folder,files(index).name), ...
        'FileType','text','Delimiter','\t');
end
summary=sortrows(vertcat(parts{:}),'taskIndex');
writetable(summary,fullfile(outputRoot,'spatial_shape_control_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

labels = replace(string(summary.control),'_',' ');
fig = figure('Color','w','Position',[100 100 1500 980]);
layout=tiledlayout(fig,2,3,'TileSpacing','loose','Padding','compact');
local_counts(nexttile(layout),labels,summary.countAbsEigLE0p05, ...
    'Near-zero eigenvalues','|\lambda| \leq 0.05');
local_counts(nexttile(layout),labels,summary.countSingularLE0p05, ...
    'Small singular values','\sigma \leq 0.05');
local_counts(nexttile(layout),labels,summary.countSchurEigLE0p05, ...
    'E-I Schur branch','|\lambda(H)| \leq 0.05');
local_counts(nexttile(layout),labels,summary.maximumRealEigenvalue, ...
    'Leading real eigenvalue','max Re(\lambda)');
local_counts(nexttile(layout),labels,summary.stableRank, ...
    'Stable rank','||J||_F^2 / ||J||_2^2');
local_counts(nexttile(layout),labels,summary.henriciDeparture, ...
    'Henrici departure from normality','normalized departure');
title(layout,'Spatial footprint controls at fixed operating point and row sums');
exportgraphics(fig,fullfile(outputRoot,'07_spatial_shape_spectral_controls.pdf'), ...
    'ContentType','vector');
savefig(fig,fullfile(outputRoot,'07_spatial_shape_spectral_controls.fig'));

selected = {'baseline','joint_left_half','joint_triangle'};
fig2=figure('Color','w','Position',[100 100 1320 940]);
layout2=tiledlayout(fig2,2,3,'TileSpacing','loose','Padding','compact');
for index=1:numel(selected)
    loaded=load(fullfile(outputRoot,['shape_' selected{index} '.mat']), ...
        'l4Kernel','l6Kernel');
    local_kernel(nexttile(layout2),loaded.l4Kernel, ...
        [replace(selected{index},'_',' ') ' - L4'],false);
    local_kernel(nexttile(layout2,index+3),loaded.l6Kernel, ...
        [replace(selected{index},'_',' ') ' - L6'],true);
end
title(layout2,'Normalized spatial footprints used by the Jacobian controls');
exportgraphics(fig2,fullfile(outputRoot,'08_spatial_shape_kernel_footprints.pdf'), ...
    'ContentType','vector');
savefig(fig2,fullfile(outputRoot,'08_spatial_shape_kernel_footprints.fig'));
end

function local_counts(ax,labels,values,panelTitle,yLabel)
bar(ax,values,'FaceColor',[0.20 0.55 0.75]);
ax.XTick=1:numel(labels); ax.XTickLabel=labels; ax.XTickLabelRotation=28;
ylabel(ax,yLabel); title(ax,panelTitle); grid(ax,'on'); box(ax,'on');
end

function local_kernel(ax,kernel,panelTitle,showXLabel)
x = -floor(size(kernel,2)/2):(ceil(size(kernel,2)/2)-1);
y = -floor(size(kernel,1)/2):(ceil(size(kernel,1)/2)-1);
imagesc(ax,x,y,kernel); axis(ax,'image'); title(ax,panelTitle);
if showXLabel; xlabel(ax,'horizontal offset'); end
ylabel(ax,'vertical offset');
colormap(ax,parula); colorbar(ax);
end
