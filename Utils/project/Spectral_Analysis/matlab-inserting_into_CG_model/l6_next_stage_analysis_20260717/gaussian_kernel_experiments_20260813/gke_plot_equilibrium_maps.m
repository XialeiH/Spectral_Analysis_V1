function gke_plot_equilibrium_maps(outputRoot)
% Plot all manipulation-specific equilibria with shared population scales.

catalog=gke_condition_catalog();
data=cell(height(catalog),1);
for index=1:height(catalog)
    file=fullfile(outputRoot,'equilibrium_only','data',catalog.name(index)+'.mat');
    if ~isfile(file); error('GKE:MissingEquilibrium','Missing %s.',file); end
    data{index}=load(file,'condition','metadata','equilibrium','summary');
end
figureRoot=fullfile(outputRoot,'equilibrium_only','figures');
if ~exist(figureRoot,'dir'); mkdir(figureRoot); end
populations=["S" "C" "I" "weighted E"];
limits=local_limits(data,populations);

figureHandle=figure('Visible','off','Color','w','Position',[50 50 1500 2500]);
layout=tiledlayout(height(catalog),numel(populations),'TileSpacing','compact', ...
    'Padding','compact');
for index=1:height(catalog)
    for populationIndex=1:numel(populations)
        nexttile;
        values=local_population(data{index}.equilibrium,populations(populationIndex), ...
            data{index}.metadata.CWeight);
        imagesc(reshape(values,40,40),limits(populationIndex,:)); axis image off;
        colormap(gca,parula); colorbar;
        if populationIndex==1
            ylabel(strrep(data{index}.condition.label,'_',' '), ...
                'Interpreter','none','FontSize',7);
        end
        if index==1; title(populations(populationIndex)); end
    end
end
title(layout,'Manipulation-specific equilibrium firing-rate maps (shared scale per population)');
local_save(figureHandle,figureRoot,'all_manipulation_equilibrium_maps');

for populationIndex=1:numel(populations)
    figureHandle=figure('Visible','off','Color','w','Position',[50 50 1500 1200]);
    layout=tiledlayout(4,5,'TileSpacing','compact','Padding','compact');
    for index=1:height(catalog)
        nexttile;
        values=local_population(data{index}.equilibrium,populations(populationIndex), ...
            data{index}.metadata.CWeight);
        imagesc(reshape(values,40,40),limits(populationIndex,:)); axis image;
        colormap(gca,parula); colorbar;
        title(data{index}.condition.label,'Interpreter','none','FontSize',8);
        xlabel('x pixel'); ylabel('y pixel');
    end
    title(layout,populations(populationIndex)+ ...
        " manipulation-specific equilibrium maps; shared color scale");
    local_save(figureHandle,figureRoot, ...
        "all_conditions_"+replace(lower(populations(populationIndex))," ","_")+ ...
        "_equilibrium_maps");
end

for index=1:height(catalog)
    figureHandle=figure('Visible','off','Color','w','Position',[100 100 1400 360]);
    layout=tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
    for populationIndex=1:numel(populations)
        nexttile;
        values=local_population(data{index}.equilibrium,populations(populationIndex), ...
            data{index}.metadata.CWeight);
        imagesc(reshape(values,40,40),limits(populationIndex,:)); axis image;
        colormap(gca,parula); colorbar; title(populations(populationIndex));
        xlabel('x pixel'); ylabel('y pixel');
    end
    title(layout,data{index}.condition.label+ ...
        sprintf(' equilibrium; residual %.2e',data{index}.summary.fixedPointResidual), ...
        'Interpreter','none');
    local_save(figureHandle,figureRoot,data{index}.condition.name+"_equilibrium_maps");
end
end

function limits=local_limits(data,populations)
limits=zeros(numel(populations),2);
for populationIndex=1:numel(populations)
    values=[];
    for index=1:numel(data)
        values=[values;local_population(data{index}.equilibrium, ...
            populations(populationIndex),data{index}.metadata.CWeight)]; %#ok<AGROW>
    end
    limits(populationIndex,:)=[min(values) max(values)];
end
end

function values=local_population(equilibrium,population,cWeight)
switch population
    case "S"; values=equilibrium.S(:);
    case "C"; values=equilibrium.C(:);
    case "I"; values=equilibrium.I(:);
    case "weighted E"; values=(1-cWeight)*equilibrium.S(:)+cWeight*equilibrium.C(:);
end
end

function local_save(figureHandle,root,name)
name=string(name);
exportgraphics(figureHandle,fullfile(root,name+".pdf"),'ContentType','vector');
exportgraphics(figureHandle,fullfile(root,name+".png"),'Resolution',200);
savefig(figureHandle,fullfile(root,name+".fig"));
close(figureHandle);
end
