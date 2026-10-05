function gke_plot_condition(dataFile,figureRoot)
% Generate report-ready figures for one completed manipulation.

loaded=load(dataFile,'condition','metadata','operators','equilibrium', ...
    'eigenvalues','singularValues','singularModes');
condition=loaded.condition;
tag=sprintf('%02d_%s',condition.id,condition.name);
conditionRoot=fullfile(figureRoot,'conditions',tag);
if ~exist(conditionRoot,'dir'); mkdir(conditionRoot); end

figureHandle=figure('Color','w','Position',[80 80 1650 480]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile
scatter(real(loaded.eigenvalues),imag(loaded.eigenvalues),7,'.');
xline(1,'r--','LineWidth',1); xline(0,'k:'); yline(0,'k:');
xlabel('Re(\lambda)'); ylabel('Im(\lambda)'); grid on
title('Complete eigenspectrum');
nexttile
edges=(floor(min(real(loaded.eigenvalues))/0.05): ...
    ceil(max(real(loaded.eigenvalues))/0.05))*0.05;
histogram(real(loaded.eigenvalues),edges,'FaceColor',[0.15 0.45 0.75], ...
    'EdgeColor','none');
xlabel('Re(\lambda)'); ylabel('Count'); grid on
title('Eigenvalue histogram; bin width 0.05');
nexttile
semilogy(1:numel(loaded.singularValues),loaded.singularValues, ...
    'Color',[0.15 0.45 0.75],'LineWidth',1.2);
xlabel('Descending singular-value index'); ylabel('\sigma_k'); grid on
title('Complete singular spectrum');
sgtitle(condition.label+sprintf(' | %s',loaded.metadata.Definition), ...
    'Interpreter','none','FontWeight','bold','FontSize',13);
local_save(figureHandle,conditionRoot,'01_spectrum_histogram_and_singular_values');

populations=["S" "C" "I" "weighted E"];
figureHandle=figure('Color','w','Position',[80 80 1580 440]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
for populationIndex=1:4
    nexttile
    values=local_population(loaded.equilibrium,populations(populationIndex), ...
        loaded.metadata.CWeight);
    imagesc(reshape(values,[40 40])); axis image off
    colormap(gca,parula); colorbar
    title(populations(populationIndex)+" firing rate");
end
sgtitle(condition.label+" at its manipulation-specific equilibrium", ...
    'Interpreter','none','FontWeight','bold','FontSize',14);
local_save(figureHandle,conditionRoot,'02_equilibrium_firing_rate_maps');

figureHandle=figure('Color','w','Position',[40 40 1580 820]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
rightVector=loaded.singularModes.LowRightVectors(:,end);
leftVector=loaded.singularModes.LowLeftVectors(:,end);
for row=1:2
    if row==1; vector=rightVector; role="Right/input"; else; vector=leftVector; role="Left/output"; end
    for populationIndex=1:4
        nexttile
        values=local_mode_population(vector,populations(populationIndex), ...
            loaded.metadata.CWeight);
        bound=max(abs(values)); if bound==0; bound=1; end
        imagesc(reshape(values,[40 40])); axis image off
        colormap(gca,turbo); clim([-bound bound]); colorbar
        title(role+" "+populations(populationIndex));
    end
end
sgtitle(sprintf('%s | smallest singular mode, \\sigma_{min}=%.3e', ...
    condition.label,loaded.singularModes.LowValues(end)), ...
    'Interpreter','tex','FontWeight','bold','FontSize',14);
local_save(figureHandle,conditionRoot,'03_smallest_left_and_right_singular_vectors');

figureHandle=figure('Color','w','Position',[80 80 1050 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for index=1:2
    nexttile
    if index==1; matrix=loaded.operators.C_SS; titleText='Representative E-source kernel: C_{SS}';
    else; matrix=loaded.operators.C_SI; titleText='Representative I-source kernel: C_{SI}'; end
    kernel=local_centered_row(matrix,[40 40]);
    imagesc(kernel); axis image
    colormap(gca,parula); colorbar
    title(titleText,'Interpreter','tex'); xlabel('x offset'); ylabel('y offset');
end
sgtitle(condition.label+" | representative L4 spatial footprints", ...
    'Interpreter','none','FontWeight','bold','FontSize',14);
local_save(figureHandle,conditionRoot,'04_representative_E_and_I_kernel_footprints');
end

function values=local_population(equilibrium,population,cWeight)
switch population
    case "S"; values=equilibrium.S(:);
    case "C"; values=equilibrium.C(:);
    case "I"; values=equilibrium.I(:);
    case "weighted E"; values=(1-cWeight)*equilibrium.S(:)+cWeight*equilibrium.C(:);
end
end

function values=local_mode_population(vector,population,cWeight)
n=numel(vector)/3;
switch population
    case "S"; values=real(vector(1:n));
    case "C"; values=real(vector(n+(1:n)));
    case "I"; values=real(vector(2*n+(1:n)));
    case "weighted E"
        values=real((1-cWeight)*vector(1:n)+cWeight*vector(n+(1:n)));
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
