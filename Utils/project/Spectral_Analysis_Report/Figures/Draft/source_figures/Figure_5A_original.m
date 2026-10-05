function Figure_5A
% Combine the compensation heatmap and two fixed-point maps.

projectRoot = [repro_paths('project') ''];
modelRoot = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model', 'Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main');
analysisRoot = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model', 'l6_next_stage_analysis_20260717');
codeRoot = fullfile(analysisRoot, 'l6_ei_compensation_grid_20260803');
runtimeRoot = fullfile(analysisRoot, 'real_tuning_l6_i_20260722', ...
    'mechanism_analysis_20260722', 'runtime_h96');
setupFile = fullfile(analysisRoot, 'results_global_bifurcation_20260722', ...
    'global_bifurcation_setup.mat');
sourceFig = fullfile(modelRoot, 'Figures', ...
    'spectral_analysis_eigenvalue_eigenvectors', 'L6 and Inhibition', ...
    'L6_and_EI_Inhibition_offset', ...
    '6.1.2_FullGrid_Beta6_BetaEI_Canonical_HCnorm_Heatmap_with_HCnorm3_Boundary.fig');
outputDir = fullfile(projectRoot, 'Spectral_Analysis_Report', ...
    'Figures', 'Final draft');
outputPdf = fullfile(outputDir, 'Figure 5A.pdf');
outputFig = fullfile(outputDir, 'Figure 5A.fig');

% Compute the same two fixed-point conditions used by Figure 5A.2.
repro_addpath(modelRoot);
repro_addpath(fullfile(modelRoot, 'Utils'));
repro_addpath(runtimeRoot, '-begin');
repro_addpath(codeRoot, '-begin');
loaded = load(setupFile, 'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
regressionSlope = 0.0729883503336975;
beta6 = [0, 0.35];
betaEI = [0, regressionSlope * beta6(2)];
states = nan(numel(baseline), 2);
hcNorm = nan(1, 2);
for condition = 1:2
    phi = @(state) l6ns_phi_l6_ei(state, 1 + beta6(condition), ...
        1 + betaEI(condition), context);
    fixed = real_tuning_fixed_point(phi, baseline, context.RelaxationP);
    assert(fixed.Converged, 'Condition %d did not converge.', condition);
    states(:, condition) = fixed.State(:);
    hcNorm(condition) = HC_norm_diff( ...
        fixed.State(1:n), fixed.State(n+(1:n)), ...
        fixed.State(2*n+(1:n)), baseline(1:n), baseline(n+(1:n)), ...
        baseline(2*n+(1:n)), context.CWeight, 0.8, 0.2);
end
eMaps = cell(1, 2);
for condition = 1:2
    sMap = reshape(states(1:n, condition), mapSize);
    cMap = reshape(states(n+(1:n), condition), mapSize);
    eMaps{condition} = (1-context.CWeight)*sMap + context.CWeight*cMap;
end
colorLimits = [min([eMaps{1}(:); eMaps{2}(:)]), ...
    max([eMaps{1}(:); eMaps{2}(:)])];

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.3 0.3 22 12], 'Renderer', 'painters');

% Copy the original compensation heatmap into the large left panel.
sourceHandle = openfig(sourceFig, 'invisible');
sourceCleanup = onCleanup(@() close(sourceHandle)); %#ok<NASGU>
sourceAxes = findall(sourceHandle, 'Type', 'axes');
assert(numel(sourceAxes) == 1, 'Expected one source heatmap axes.');
sourceAxes = sourceAxes(1);
leftAxes = axes(fig, 'Position', [0.070 0.115 0.430 0.788]);
copyobj(flipud(allchild(sourceAxes)), leftAxes);
set(leftAxes, 'XLim', sourceAxes.XLim, 'YLim', sourceAxes.YLim, ...
    'CLim', sourceAxes.CLim, 'ColorScale', sourceAxes.ColorScale, ...
    'YDir', sourceAxes.YDir, 'FontName', 'Arial', 'FontSize', 22, ...
    'FontWeight', 'bold', 'LineWidth', 1.5, 'TickDir', 'out', ...
    'Box', 'on', 'Layer', 'bottom');
pbaspect(leftAxes, [1 1 1]);
colormap(leftAxes, jet(256));
xlabel(leftAxes, sourceAxes.XLabel.String, ...
    'Interpreter', sourceAxes.XLabel.Interpreter, 'FontSize', 24, ...
    'FontWeight', 'bold');
ylabel(leftAxes, sourceAxes.YLabel.String, ...
    'Interpreter', sourceAxes.YLabel.Interpreter, 'FontSize', 24, ...
    'FontWeight', 'bold');
leftAxes.XLabel.Units = 'normalized';
leftAxes.XLabel.Position = [0.50 -0.070 0];
leftAxes.YLabel.Units = 'normalized';
leftAxes.YLabel.Position = [-0.072 0.50 0];
leftAxes.Toolbar.Visible = 'off';

% Update displayed units and add the two selected-condition markers.
displayObjects = findall(leftAxes);
for index = 1:numel(displayObjects)
    if isprop(displayObjects(index), 'DisplayName')
        displayObjects(index).DisplayName = replace( ...
            string(displayObjects(index).DisplayName), 'Hz', 'sp/s');
    end
end
regressionLine = findobj(leftAxes, '-regexp', 'DisplayName', ...
    '^minimum-HC regression');
assert(~isempty(regressionLine), 'Minimum-HC regression line not found.');
regressionLine = regressionLine(1);
xSecond = 0.35;
ySecond = interp1(regressionLine.XData, regressionLine.YData, ...
    xSecond, 'linear');
hold(leftAxes, 'on');
markerWidth = 0.0065;
markerHeight = 0.0050;
originMarker = patch(leftAxes, [0 markerWidth markerWidth 0], ...
    [0 0 markerHeight markerHeight], 'w', 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');
secondMarker = patch(leftAxes, ...
    [xSecond-markerWidth xSecond xSecond xSecond-markerWidth], ...
    ySecond + 0.5*markerHeight*[-1 -1 1 1], 'w', ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
uistack([originMarker secondMarker], 'top');

legendObjects = gobjects(1, 5);
patterns = {'^moved-FP HC slope', '^minimum-HC regression', ...
    '^HC norm = 2', '^HC norm = 5', '^HC norm > 170'};
for index = 1:numel(patterns)
    object = findobj(leftAxes, '-regexp', 'DisplayName', patterns{index});
    assert(~isempty(object), 'Missing legend object %s.', patterns{index});
    legendObjects(index) = object(1);
end
legendHandle = legend(leftAxes, legendObjects, ...
    get(legendObjects, {'DisplayName'}), 'NumColumns', 3, ...
    'Box', 'off', 'FontName', 'Arial', 'FontSize', 19, ...
    'FontWeight', 'bold');
legendHandle.Units = 'normalized';
legendHandle.Position = [0.070 0.918 0.485 0.065];

leftColorbar = colorbar(leftAxes, 'eastoutside');
leftColorbar.Units = 'normalized';
leftColorbar.Position = [0.515 0.115 0.015 0.788];
leftColorbar.FontName = 'Arial';
leftColorbar.FontSize = 20;
leftColorbar.FontWeight = 'bold';
leftColorbar.Title.String = 'sp/s';
leftColorbar.Title.FontSize = 22;
leftColorbar.Title.FontWeight = 'bold';
leftColorbar.Label.String = 'HC norm';
leftColorbar.Label.FontSize = 22;
leftColorbar.Label.FontWeight = 'bold';

% Stack the two fixed-point maps on the right.
rightPositions = [0.675 0.582 0.175 0.321; ...
                  0.675 0.115 0.175 0.321];
rightAxes = gobjects(1, 2);
for condition = 1:2
    ax = axes(fig, 'Position', rightPositions(condition, :));
    rightAxes(condition) = ax;
    imagesc(ax, eMaps{condition});
    axis(ax, 'image');
    set(ax, 'YDir', 'normal', 'FontName', 'Arial', 'FontSize', 20, ...
        'FontWeight', 'bold', 'LineWidth', 1.4, 'TickDir', 'out');
    clim(ax, colorLimits);
    colormap(ax, jet(256));
    xticks(ax, [1 20 40]);
    yticks(ax, [1 20 40]);
    xlabel(ax, 'Horizontal pixel index (pixel)', ...
        'FontSize', 21, 'FontWeight', 'bold');
    ylabel(ax, 'Vertical pixel index (pixel)', ...
        'FontSize', 21, 'FontWeight', 'bold');
    if condition == 1
        panelTitle = 'Condition 1';
    else
        panelTitle = {'Condition 2', ...
            sprintf('HCnorm = %.3f sp/s', hcNorm(condition))};
    end
    title(ax, panelTitle, 'FontSize', 22, 'FontWeight', 'bold');
    ax.Toolbar.Visible = 'off';
end

rightColorbar = colorbar(rightAxes(1), 'eastoutside');
rightColorbar.Units = 'normalized';
rightColorbar.Position = [0.865 0.582 0.014 0.321];
rightColorbar.FontName = 'Arial';
rightColorbar.FontSize = 20;
rightColorbar.FontWeight = 'bold';
rightColorbar.Title.String = 'sp/s';
rightColorbar.Title.FontSize = 22;
rightColorbar.Title.FontWeight = 'bold';
rightColorbar.Label.String = 'E firing rate';
rightColorbar.Label.FontSize = 22;
rightColorbar.Label.FontWeight = 'bold';

savefig(fig, outputFig);
if isfile(outputPdf)
    delete(outputPdf);
end
exportgraphics(fig, outputPdf, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
close(fig);
fprintf('Saved %s and %s.\n', outputPdf, outputFig);
end
