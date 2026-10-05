function plot_l6ns_branch_atlas_phase1(root)
% Plot reconciled branch geometry, stability, and singularity diagnostics.

if nargin<1 || isempty(root); root=pwd; end
loaded=load(fullfile(root,'branch_atlas_20260803','aggregate', ...
    'phase1_branch_atlas.mat'),'atlas');
atlas=loaded.atlas;
outputDir=fullfile(root,'branch_atlas_20260803','figures');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
pathways=["L6","Inhibition"];
for pathwayIndex=1:numel(pathways)
    pathway=pathways(pathwayIndex);
    points=atlas.Points(atlas.Points.pathway==pathway,:);
    figureHandle=figure('Color','w','Position',[80 80 1500 900]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    nexttile;
    scatter(points.freezeWeight,points.distanceFromPersistent,18, ...
        points.unstableDimension,'filled');
    xlabel('freeze weight w'); ylabel('||f-f_*||_2');
    title(sprintf('%s equilibrium components',pathway));
    colorbar; grid on; box on;
    nexttile;
    scatter(points.freezeWeight,points.maxRealLambda,18, ...
        points.unstableDimension,'filled'); hold on;
    yline(1,'r--','stability boundary');
    xlabel('freeze weight w'); ylabel('max Re lambda(J)');
    title('Spectral stability and Morse index'); colorbar; grid on; box on;
    nexttile;
    semilogy(points.freezeWeight,max(points.sigmaMinA,eps),'.','MarkerSize',8); hold on;
    semilogy(points.freezeWeight,max(points.sigma2A,eps),'.','MarkerSize',8);
    xlabel('freeze weight w'); ylabel('singular value of J-I');
    legend({'sigma_{min}','sigma_2'},'Location','best');
    title('Steady-state singularity diagnostics'); grid on; box on;
    nexttile;
    scatter(points.freezeWeight,points.dwds,18,points.componentId,'filled'); hold on;
    yline(0,'k--');
    xlabel('freeze weight w'); ylabel('dw/ds');
    title('Arclength turning geometry'); colorbar; grid on; box on;
    sgtitle(sprintf('%s Phase I branch atlas',pathway));
    stem=fullfile(outputDir,sprintf('%s_phase1_branch_atlas',lower(pathway)));
    exportgraphics(figureHandle,[stem '.pdf'],'ContentType','vector');
    savefig(figureHandle,[stem '.fig']);
end

figureHandle=figure('Color','w','Position',[100 100 1300 620]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for pathwayIndex=1:numel(pathways)
    nexttile;
    pathway=pathways(pathwayIndex);
    points=atlas.Points(atlas.Points.pathway==pathway,:);
    scatter(points.meanC,points.meanI,20,points.freezeWeight,'filled');
    xlabel('mean C rate'); ylabel('mean I rate');
    title(sprintf('%s component separation',pathway));
    colorbar; grid on; box on;
end
sgtitle('Equilibrium families in population-mean coordinates');
stem=fullfile(outputDir,'phase1_component_population_projection');
exportgraphics(figureHandle,[stem '.pdf'],'ContentType','vector');
savefig(figureHandle,[stem '.fig']);
end
