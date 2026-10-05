function plot_three_hcnorm3_hcoptimal_slopes()
% Plot the three HCnorm heatmaps with HC-optimal fixed-point slopes.

figureBase = fullfile( ...
    [repro_paths('project') '/Spectral_Analysis'], ...
    'matlab-inserting_into_CG_model/Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors', ...
    'L6 and Inhibition');
slopes = audit_hc_optimal_fixedpoint_slopes([],figureBase);

local_plot_l6_inhibition(fullfile(figureBase,'L6_and_Inhibition_offset'), ...
    slopes.L6Inhibition);
local_plot_pair(fullfile(figureBase,'L6_and_Excitation_offset'), ...
    '\beta_6 (fractional whole-L6 increase)', ...
    '\beta_E (fractional recurrent-excitation change)', ...
    '08_FullGrid_Beta6_BetaE_Canonical_HCnorm_Heatmap_with_HCOptimalSlope', ...
    slopes.L6Excitation);
local_plot_pair(fullfile(figureBase,'Excitation_and_Inhibition_offset'), ...
    '\beta_E (fractional recurrent-excitation increase)', ...
    '\beta_I (fractional whole-inhibition increase)', ...
    '09_FullGrid_BetaE_BetaI_Canonical_HCnorm_Heatmap_with_HCOptimalSlope', ...
    slopes.ExcitationInhibition);
end

function local_plot_l6_inhibition(outputRoot,theorySlope)
data = load(fullfile(outputRoot,'full_grid_results.mat'), ...
    'beta6Values','betaIValues','validHC');
fig = local_heatmap(data.beta6Values,data.betaIValues,data.validHC, ...
    '\beta_6 (fractional whole-L6 increase)', ...
    '\beta_I (fractional whole-inhibition increase)');
ax = gca;
hold(ax,'on');
hTheory = local_add_theory_line(ax,data.beta6Values,data.betaIValues, ...
    theorySlope);
hBoundary = local_add_boundary(ax,data.beta6Values,data.betaIValues, ...
    data.validHC);
local_legend(ax,hTheory,hBoundary,theorySlope);
local_save(fig,outputRoot, ...
    '07_FullGrid_Beta6_BetaI_Canonical_HCnorm_Heatmap_with_HCOptimalSlope');
end

function local_plot_pair(outputRoot,xLabel,yLabel,stem,theorySlope)
data = load(fullfile(outputRoot,'full_grid_results.mat'));
if isfield(data,'betaXValues')
    xValues = data.betaXValues;
    yValues = data.betaYValues;
else
    xValues = data.betaValues;
    yValues = data.betaValues;
end
fig = local_heatmap(xValues,yValues,data.plotHC,xLabel,yLabel);
ax = gca;
hold(ax,'on');
hTheory = local_add_theory_line(ax,xValues,yValues,theorySlope);
hBoundary = local_add_boundary(ax,xValues,yValues,data.validHC);
local_legend(ax,hTheory,hBoundary,theorySlope);
local_save(fig,outputRoot,stem);
end

function local_legend(ax,hTheory,hBoundary,theorySlope)
legend(ax,[hTheory hBoundary], ...
    {sprintf('HC-optimal fixed-point slope %.6f',theorySlope), ...
    'HC norm = 3 Hz boundary'}, ...
    'Location','northoutside','Orientation','horizontal', ...
    'FontSize',15,'FontWeight','bold','Box','off');
end

function hTheory = local_add_theory_line(ax,xValues,yValues,slope)
lineX = linspace(min(xValues),max(xValues),500);
lineY = slope*lineX;
inside = lineY>=min(yValues) & lineY<=max(yValues);
hTheory = plot(ax,lineX(inside),lineY(inside),'k--', ...
    'LineWidth',2.2,'DisplayName', ...
    sprintf('HC-optimal fixed-point slope %.6f',slope));
end

function fig = local_heatmap(xValues,yValues,plotHC,xLabel,yLabel)
fig = figure('Visible','off','Color','w','Position',[80 80 1280 1020]);
heatmapImage = imagesc(xValues,yValues,plotHC);
heatmapImage.AlphaData = isfinite(plotHC);
ax = gca;
set(ax,'YDir','normal','Color',[0.80 0.80 0.80]);
axis(ax,'xy');
axis(ax,'tight');
xlabel(ax,xLabel,'FontSize',22,'FontWeight','bold');
ylabel(ax,yLabel,'FontSize',22,'FontWeight','bold');
set(ax,'FontSize',19,'FontWeight','bold','LineWidth',1.25, ...
    'XGrid','off','YGrid','off','ColorScale','log');
colormap(fig,jet);
finiteValues = sort(plotHC(isfinite(plotHC) & plotHC>0));
lowerIndex = max(1,round(0.001*(numel(finiteValues)-1))+1);
upperIndex = min(numel(finiteValues), ...
    round(0.98*(numel(finiteValues)-1))+1);
if finiteValues(lowerIndex)<finiteValues(upperIndex)
    clim(ax,[finiteValues(lowerIndex) finiteValues(upperIndex)]);
end
cb = colorbar(ax);
cb.FontSize = 17;
cb.Label.String = 'HC norm';
cb.Label.FontSize = 19;
cb.Label.FontWeight = 'bold';
cb.Title.String = 'Hz';
cb.Title.FontSize = 17;
cb.Title.FontWeight = 'bold';
end

function hBoundary = local_add_boundary(ax,xValues,yValues,validHC)
[~,hBoundary] = contour(ax,xValues,yValues,validHC,[3 3], ...
    'LineColor','k','LineWidth',1.25, ...
    'DisplayName','HC norm = 3 Hz boundary');
end

function local_save(fig,outputRoot,stem)
savefig(fig,fullfile(outputRoot,[stem '.fig']));
exportgraphics(fig,fullfile(outputRoot,[stem '.pdf']), ...
    'ContentType','vector');
close(fig);
end
