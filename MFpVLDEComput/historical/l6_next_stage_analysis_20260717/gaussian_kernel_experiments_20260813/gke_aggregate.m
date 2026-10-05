function gke_aggregate(outputRoot)
% Aggregate all manipulation-specific spectra and produce comparison figures.

if nargin<1 || isempty(outputRoot); outputRoot=getenv('GKE_OUTPUT'); end
catalog=gke_condition_catalog();
figureRoot=fullfile(outputRoot,'figures');
if ~exist(figureRoot,'dir'); mkdir(figureRoot); end
rows=cell(height(catalog),1);
data=cell(height(catalog),1);
for index=1:height(catalog)
    summaryFile=fullfile(outputRoot,'summaries',catalog.name(index)+'.tsv');
    dataFile=fullfile(outputRoot,'data',catalog.name(index)+'.mat');
    if ~isfile(summaryFile) || ~isfile(dataFile)
        error('GKE:MissingResult','Missing completed result for %s.',catalog.name(index));
    end
    rows{index}=readtable(summaryFile,'FileType','text','Delimiter','\t');
    data{index}=load(dataFile,'condition','metadata','operatorAudit', ...
        'equilibrium','operators','eigenvalues','singularValues','singularModes');
    gke_plot_condition(dataFile,figureRoot);
end
summary=vertcat(rows{:});
writetable(summary,fullfile(outputRoot,'gaussian_kernel_experiment_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
local_write_operator_audit(outputRoot,catalog,data);
local_plot_counts(figureRoot,summary);
local_plot_stability(figureRoot,summary);
local_plot_group_spectra(figureRoot,catalog,data);
local_plot_group_histograms(figureRoot,catalog,data);
local_plot_group_singular_spectra(figureRoot,catalog,data);
local_plot_equilibrium_overview(figureRoot,catalog,data);
local_plot_kernel_overview(figureRoot,catalog,data);
end

function local_write_operator_audit(outputRoot,catalog,data)
rows=cell(height(catalog),1);
for index=1:height(catalog)
    current=data{index}.operatorAudit;
    current.conditionId(:)=catalog.id(index);
    current.conditionName(:)=catalog.name(index);
    rows{index}=current;
end
audit=vertcat(rows{:});
writetable(audit,fullfile(outputRoot,'all_operator_invariant_audits.tsv'), ...
    'FileType','text','Delimiter','\t');
end

function local_plot_counts(figureRoot,summary)
figureHandle=figure('Color','w','Position',[40 40 1750 900]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile
bar([summary.eigLE0p01 summary.eigLE0p025 summary.eigLE0p05],'grouped');
ylabel('Eigenvalue count'); grid on
legend({'|\lambda| <= 0.01','|\lambda| <= 0.025','|\lambda| <= 0.05'}, ...
    'Location','northoutside','Orientation','horizontal');
title('Near-zero eigenvalue cluster after re-equilibrating every model');
nexttile
bar([summary.sigmaLE0p01 summary.sigmaLE0p025 summary.sigmaLE0p05],'grouped');
ylabel('Singular-value count'); grid on
legend({'\sigma <= 0.01','\sigma <= 0.025','\sigma <= 0.05'}, ...
    'Location','northoutside','Orientation','horizontal');
title('Small-singular-value sector');
xticks(1:height(summary)); xticklabels(summary.name); xtickangle(35);
sgtitle('Gaussian spatial-kernel controls: near-zero comparison', ...
    'FontWeight','bold','FontSize',16);
local_save(figureHandle,figureRoot,'01_near_zero_eigenvalue_and_singular_value_counts');
end

function local_plot_stability(figureRoot,summary)
figureHandle=figure('Color','w','Position',[60 60 1650 560]);
yyaxis left
bar(summary.maximumRealEigenvalue,'FaceColor',[0.15 0.45 0.75]);
ylabel('max Re(\lambda)'); yline(1,'r--','LineWidth',1.3);
yyaxis right
semilogy(1:height(summary),summary.fixedPointResidual,'ko-','LineWidth',1);
ylabel('Fixed-point residual');
xticks(1:height(summary)); xticklabels(summary.name); xtickangle(35); grid on
title('Stability and fixed-point verification for every manipulation');
local_save(figureHandle,figureRoot,'02_stability_and_fixed_point_residuals');
end

function local_plot_group_spectra(figureRoot,catalog,data)
groups=unique(catalog.group,'stable');
for groupIndex=1:numel(groups)
    indices=find(catalog.group==groups(groupIndex));
    realLimits=[Inf -Inf]; imaginaryLimit=0;
    for index=indices'
        values=data{index}.eigenvalues;
        realLimits=[min(realLimits(1),min(real(values))) max(realLimits(2),max(real(values)))];
        imaginaryLimit=max(imaginaryLimit,max(abs(imag(values))));
    end
    [rows,columns]=local_panel_shape(numel(indices));
    figureHandle=figure('Color','w','Position',[30 30 1600 420*rows]);
    tiledlayout(rows,columns,'TileSpacing','compact','Padding','compact');
    for index=indices'
        nexttile
        values=data{index}.eigenvalues;
        scatter(real(values),imag(values),6,'.'); xline(1,'r--'); yline(0,'k:');
        xlim(realLimits+[-0.03 0.03]*max(diff(realLimits),1));
        ylim([-1 1]*max(imaginaryLimit*1.05,0.1)); grid on
        xlabel('Re(\lambda)'); ylabel('Im(\lambda)'); title(catalog.label(index));
    end
    sgtitle("Complete eigenspectra: "+groups(groupIndex)+" controls", ...
        'FontWeight','bold','FontSize',15);
    local_save(figureHandle,figureRoot,sprintf('03_%s_complete_eigenspectra',groups(groupIndex)));
end
end

function local_plot_group_histograms(figureRoot,catalog,data)
groups=unique(catalog.group,'stable');
for groupIndex=1:numel(groups)
    indices=find(catalog.group==groups(groupIndex));
    [rows,columns]=local_panel_shape(numel(indices));
    figureHandle=figure('Color','w','Position',[30 30 1600 390*rows]);
    tiledlayout(rows,columns,'TileSpacing','compact','Padding','compact');
    for index=indices'
        nexttile
        values=real(data{index}.eigenvalues);
        edges=(floor(min(values)/0.05):ceil(max(values)/0.05))*0.05;
        histogram(values,edges,'EdgeColor','none'); grid on
        xlabel('Re(\lambda)'); ylabel('Count');
        title(catalog.label(index)+" | bin 0.05");
    end
    sgtitle("Eigenvalue histograms: "+groups(groupIndex)+" controls", ...
        'FontWeight','bold','FontSize',15);
    local_save(figureHandle,figureRoot,sprintf('04_%s_eigenvalue_histograms_bin0p05',groups(groupIndex)));
end
end

function local_plot_group_singular_spectra(figureRoot,catalog,data)
groups=unique(catalog.group,'stable');
for groupIndex=1:numel(groups)
    indices=find(catalog.group==groups(groupIndex));
    figureHandle=figure('Color','w','Position',[60 60 1250 650]); hold on
    for index=indices'
        values=data{index}.singularValues;
        semilogy(1:numel(values),values,'LineWidth',1.1, ...
            'DisplayName',catalog.label(index));
    end
    xlabel('Descending singular-value index'); ylabel('\sigma_k'); grid on
    legend('Location','eastoutside','Interpreter','none');
    title("Complete singular spectra: "+groups(groupIndex)+" controls", ...
        'FontWeight','bold');
    local_save(figureHandle,figureRoot,sprintf('05_%s_complete_singular_spectra',groups(groupIndex)));
end
end

function local_plot_equilibrium_overview(figureRoot,catalog,data)
figureHandle=figure('Color','w','Position',[20 20 1800 4100]);
tiledlayout(height(catalog),4,'TileSpacing','compact','Padding','compact');
populations=["S" "C" "I" "weighted E"];
limits=local_population_limits(data,populations);
for index=1:height(catalog)
    for populationIndex=1:4
        nexttile
        values=local_population(data{index}.equilibrium,populations(populationIndex), ...
            data{index}.metadata.CWeight);
        imagesc(reshape(values,[40 40])); axis image off
        colormap(gca,parula); clim(limits(populationIndex,:)); colorbar
        title(catalog.label(index)+": "+populations(populationIndex), ...
            'Interpreter','none','FontSize',7);
    end
end
sgtitle('Manipulation-specific equilibrium firing-rate maps; shared population scales', ...
    'FontWeight','bold','FontSize',16);
local_save(figureHandle,figureRoot,'06_all_equilibrium_firing_rate_maps');
end

function local_plot_kernel_overview(figureRoot,catalog,data)
figureHandle=figure('Color','w','Position',[20 20 1200 4100]);
tiledlayout(height(catalog),2,'TileSpacing','compact','Padding','compact');
for index=1:height(catalog)
    matrices={data{index}.operators.C_SS,data{index}.operators.C_SI};
    labels=["E-source C_{SS}" "I-source C_{SI}"];
    for column=1:2
        nexttile
        imagesc(local_centered_row(matrices{column},[40 40])); axis image off
        colormap(gca,parula); colorbar
        title(catalog.label(index)+": "+labels(column), ...
            'Interpreter','tex','FontSize',7);
    end
end
sgtitle('Representative controlled L4 spatial footprints', ...
    'FontWeight','bold','FontSize',16);
local_save(figureHandle,figureRoot,'07_all_representative_E_and_I_kernel_footprints');
end

function [rows,columns]=local_panel_shape(count)
columns=min(3,count); rows=ceil(count/columns);
end

function limits=local_population_limits(data,populations)
limits=[Inf(4,1) -Inf(4,1)];
for index=1:numel(data)
    for populationIndex=1:4
        values=local_population(data{index}.equilibrium,populations(populationIndex), ...
            data{index}.metadata.CWeight);
        limits(populationIndex,:)=[min(limits(populationIndex,1),min(values)) ...
            max(limits(populationIndex,2),max(values))];
    end
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

function kernel=local_centered_row(matrix,mapSize)
kernel=reshape(full(matrix(1,:)),mapSize);
kernel=circshift(kernel,[floor(mapSize(1)/2) floor(mapSize(2)/2)]);
kernel=kernel/max(sum(abs(kernel(:))),eps);
end

function local_save(figureHandle,root,name)
exportgraphics(figureHandle,fullfile(root,[name '.pdf']),'ContentType','vector');
exportgraphics(figureHandle,fullfile(root,[name '.png']),'Resolution',180);
savefig(figureHandle,fullfile(root,[name '.fig']));
close(figureHandle)
end
