function make_Final_Figure_5E
% Final Figure 5E: pathway maps, tuning curves, and return trajectory.

projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project';
artifactRoot = fileparts(mfilename('fullpath'));
dataFile = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model', 'Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main', 'Figures', ...
    'spectral_analysis_eigenvalue_eigenvectors', 'L6 and Inhibition', ...
    'L6_and_Inhibition_offset', 'Figure6_0_FourContrast_FiveMetric', ...
    'figure6_0_contrast100_dataset.mat');
correctedTuningFile = fullfile(artifactRoot, ...
    'Figure_5E_2_corrected_tuning.mat');
trajectoryFile = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model', 'report_figure_artifacts', ...
    '1_Neural_Network_Surrogate', 'Figure_1I', ...
    'Figure_1I_trajectories.mat');
outputRoot = fullfile(projectRoot, 'Spectral_Analysis_Report', ...
    'Figures', 'Final draft');
outputPdf = fullfile(outputRoot, 'Figure 5E.pdf');
outputFig = fullfile(outputRoot, 'Figure 5E.fig');

if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

loaded = load(dataFile, 'betaGrid', 'fixedPointStates', ...
    'leadingEigenclusterEnvelopes', 'leadingRightSingularModes', ...
    'setup');
correctedTuning = load(correctedTuningFile, 'tuningCurves', ...
    'fullAngles');
trajectory = load(trajectoryFile, 'timesMs', 'trajectories');

betaGrid = double(loaded.betaGrid(:).');
targetBeta6 = [0, 0.075, 0.30];
targetBetaI = 0.561654 * targetBeta6;
pointIndices = zeros(1, 3);
for conditionIndex = 2:3
    [~, beta6Index] = min(abs(betaGrid - targetBeta6(conditionIndex)));
    [~, betaIIndex] = min(abs(betaGrid - targetBetaI(conditionIndex)));
    pointIndices(conditionIndex) = betaIIndex + ...
        (beta6Index - 1) * numel(betaGrid);
end

mapSide = 40;
populationSize = mapSide^2;
cWeight = double(loaded.setup.CWeight);
fixedMaps = zeros(mapSide, mapSide, 3);
clusterMaps = zeros(mapSide, mapSide, 3);
singularMaps = zeros(mapSide, mapSide, 3);

fixedMaps(:, :, 1) = local_e_map(loaded.setup.BaselineState, ...
    mapSide, populationSize, cWeight);
clusterMaps(:, :, 1) = local_e_map(loaded.setup.BaselineClusterEnvelope, ...
    mapSide, populationSize, cWeight);
singularMaps(:, :, 1) = abs(local_e_map(loaded.setup.BaselineRightMode, ...
    mapSide, populationSize, cWeight));

for conditionIndex = 2:3
    pointIndex = pointIndices(conditionIndex);
    fixedMaps(:, :, conditionIndex) = local_e_map( ...
        loaded.fixedPointStates(:, pointIndex), ...
        mapSide, populationSize, cWeight);
    clusterMaps(:, :, conditionIndex) = local_e_map( ...
        loaded.leadingEigenclusterEnvelopes(:, pointIndex), ...
        mapSide, populationSize, cWeight);
    singularVector = double(loaded.leadingRightSingularModes(:, pointIndex));
    if dot(singularVector, double(loaded.setup.BaselineRightMode)) < 0
        singularVector = -singularVector;
    end
    singularMaps(:, :, conditionIndex) = abs(local_e_map( ...
        singularVector, mapSide, populationSize, cWeight));
end

fixedLimits = local_limits(fixedMaps);
clusterLimits = local_limits(clusterMaps);
singularLimits = [0, max(singularMaps, [], 'all')];
tuningCurves = double(correctedTuning.tuningCurves);
tuningLimits = [floor(min(tuningCurves, [], 'all')) - 1, ...
    ceil(max(tuningCurves, [], 'all')) + 1];

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172
] / 255;
pixelColors = neuroColors([6, 1], :);
trajectoryColors = neuroColors([6, 1, 4], :);

fig = figure('Visible', 'on', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.2, 0.2, 36, 24], 'Renderer', 'painters');
columnLeft = [0.060, 0.270, 0.480, 0.705];
columnWidth = 0.150;
rowBottom = [0.690, 0.385, 0.080];
rowHeight = 0.225;
rowLabels = {'Condition 1 (Baseline)', 'Condition 2', 'Condition 3'};
columnTitles = {'E Fixed Point', 'Top Eigencluster', 'Top Singular Mode'};
mapCollections = {fixedMaps, clusterMaps, singularMaps};
mapLimits = {fixedLimits, clusterLimits, singularLimits};
heatAxes = gobjects(3, 3);

for rowIndex = 1:3
    for columnIndex = 1:3
        ax = axes(fig, 'Position', [columnLeft(columnIndex), ...
            rowBottom(rowIndex), columnWidth, rowHeight]);
        heatAxes(rowIndex, columnIndex) = ax;
        imagesc(ax, mapCollections{columnIndex}(:, :, rowIndex));
        axis(ax, 'image');
        set(ax, 'YDir', 'normal', 'FontName', 'Arial', 'FontSize', 26, ...
            'LineWidth', 1.4, 'Box', 'on', 'TickDir', 'out', ...
            'XTick', [1 10 20 30 40], 'YTick', [1 10 20 30 40]);
        clim(ax, mapLimits{columnIndex});
        colormap(ax, jet(256));
        ax.Toolbar.Visible = 'off';
        if rowIndex == 1
            title(ax, columnTitles{columnIndex}, ...
                'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
        end
        if rowIndex == 3
            xlabel(ax, 'Map column (pixel)', 'FontSize', 26);
        else
            ax.XTickLabel = [];
        end
        if columnIndex == 1
            ylabel(ax, rowLabels{rowIndex}, 'FontSize', 28, ...
                'FontWeight', 'bold');
        else
            ax.YTickLabel = [];
        end
    end
end

for columnIndex = 1:3
    ax = heatAxes(1, columnIndex);
    colorbarHandle = colorbar(ax);
    ax.Position = [columnLeft(columnIndex), rowBottom(1), ...
        columnWidth, rowHeight];
    colorbarHandle.Units = 'normalized';
    colorbarHandle.Position = [columnLeft(columnIndex) + ...
        columnWidth + 0.007, rowBottom(1), 0.011, rowHeight];
    colorbarHandle.FontName = 'Arial';
    colorbarHandle.FontSize = 26;
    colorbarHandle.LineWidth = 1.2;
    colorbarHandle.Title.String = 'sp/s';
    colorbarHandle.Title.FontSize = 26;
    colorbarHandle.Title.FontWeight = 'bold';
end

% Row 2, column 4: preserve the original six tuning curves.
angles = double(correctedTuning.fullAngles(:));
tuningAxes = axes(fig, 'Position', [columnLeft(4), rowBottom(2), ...
    columnWidth, rowHeight]);
hold(tuningAxes, 'on');
tuningLines = gobjects(6, 1);
lineIndex = 0;
lineStyles = {'--', '-', '-.'};
for pixelIndex = 1:2
    for conditionIndex = 1:3
        lineIndex = lineIndex + 1;
        tuningLines(lineIndex) = plot(tuningAxes, angles, ...
            tuningCurves(:, pixelIndex, conditionIndex), ...
            lineStyles{conditionIndex}, ...
            'Color', pixelColors(pixelIndex, :), ...
            'LineWidth', 3.0 + 0.2 * (conditionIndex > 1));
    end
end
set(tuningAxes, 'FontName', 'Arial', 'FontSize', 26, ...
    'LineWidth', 1.4, 'Box', 'on', 'TickDir', 'out', ...
    'XLim', [0 180], 'XTick', 0:45:180, 'YLim', tuningLimits);
tuningAxes.Toolbar.Visible = 'off';
tuningAxes.YGrid = 'on';
tuningAxes.XGrid = 'off';
tuningAxes.GridAlpha = 0.16;
ylabel(tuningAxes, 'E firing rate (sp/s)', 'FontSize', 26);
xlabel(tuningAxes, 'Orientation (deg)', 'FontSize', 26);

legendLabels = {'(5,10), 0 deg: baseline', ...
    '(5,10), 0 deg: Condition 2', '(5,10), 0 deg: Condition 3', ...
    '(1,10), 22.5 deg: baseline', ...
    '(1,10), 22.5 deg: Condition 2', ...
    '(1,10), 22.5 deg: Condition 3'};
legendHandle = legend(tuningAxes, tuningLines, legendLabels, ...
    'NumColumns', 2, 'FontName', 'Arial', 'FontSize', 13, ...
    'Box', 'off', 'Location', 'northoutside');
legendHandle.Units = 'normalized';
legendHandle.Position(1) = columnLeft(4) + ...
    (columnWidth - legendHandle.Position(3)) / 2;
legendHandle.Position(2) = rowBottom(2) + rowHeight + 0.010;

% Row 3, column 4: reproduce Figure 1I at the same size as a heatmap.
trajectoryAxes = axes(fig, 'Position', [columnLeft(4), rowBottom(3), ...
    columnWidth, rowHeight]);
hold(trajectoryAxes, 'on');
trajectoryLines = gobjects(3, 1);
for conditionIndex = 1:3
    trajectoryLines(conditionIndex) = plot(trajectoryAxes, ...
        trajectory.timesMs, trajectory.trajectories(:, conditionIndex), ...
        '-', 'Color', trajectoryColors(conditionIndex, :), ...
        'LineWidth', 3.2);
end
xlim(trajectoryAxes, [0 180]);
ylim(trajectoryAxes, [0, 1.06 * max(trajectory.trajectories, [], 'all')]);
xlabel(trajectoryAxes, 'Time (ms)', 'FontSize', 26);
ylabel(trajectoryAxes, 'HC norm (sp/s)', 'FontSize', 26);
set(trajectoryAxes, 'FontName', 'Arial', 'FontSize', 26, ...
    'LineWidth', 1.4, 'TickDir', 'out', 'Box', 'on', ...
    'XTick', 0:30:180, 'Layer', 'top');
trajectoryAxes.Toolbar.Visible = 'off';
grid(trajectoryAxes, 'on');
trajectoryAxes.GridColor = [0.82 0.82 0.82];
trajectoryAxes.GridAlpha = 0.55;
legend(trajectoryAxes, trajectoryLines, ...
    {'Condition 1 (Baseline)', 'Condition 2', 'Condition 3'}, ...
    'Location', 'northeast', 'FontName', 'Arial', ...
    'FontSize', 15, 'Box', 'off');

if isfile(outputPdf)
    delete(outputPdf);
end
if isfile(outputFig)
    delete(outputFig);
end
exportgraphics(fig, outputPdf, 'ContentType', 'vector', ...
    'BackgroundColor', 'white');
savefig(fig, outputFig, 'compact');
fprintf('Saved %s\nSaved %s\n', outputPdf, outputFig);
end

function map = local_e_map(vector, mapSide, populationSize, cWeight)
vector = double(vector(:));
eVector = (1 - cWeight) * vector(1:populationSize) + ...
    cWeight * vector(populationSize + (1:populationSize));
map = reshape(eVector, mapSide, mapSide);
end

function limits = local_limits(values)
finiteValues = values(isfinite(values));
limits = [min(finiteValues), max(finiteValues)];
if limits(1) == limits(2)
    limits = limits + [-1, 1] * max(1, abs(limits(1))) * 1e-6;
end
end
