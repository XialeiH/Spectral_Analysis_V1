function generate_report_weight_sweep_histograms()
% Generate matched L6- and inhibition-weight eigenvalue histogram panels.

% Each panel uses bins of width 0.05 and an x-axis adapted to that panel's
% finite real-eigenvalue range.

sourceFile = [[repro_paths('project') '/'] ...
    'Spectral_Analysis/matlab-inserting_into_CG_model/' ...
    'l6_next_stage_analysis_20260717/near_zero_clustering_20260722/' ...
    'results_canonical/near_zero_canonical_result.mat'];
l6FigureRoot = [[repro_paths('project') '/'] ...
    'Spectral_Analysis/matlab-inserting_into_CG_model/' ...
    'Complete_Code_for_Paper3/NYU-Vision-2Drive-main/Figures/' ...
    'spectral_analysis_eigenvalue_eigenvectors/eigenspectrum_L6weight'];
outputRoot = [repro_paths('project') '/Spectral_Analysis_Report'];
binWidth = 0.05;

if ~isfile(sourceFile)
    error('Missing canonical sweep file: %s', sourceFile);
end
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

loaded = load(sourceFile, 'result');
weightTable = loaded.result.WeightSweep;
spectra = loaded.result.SweepSpectra;

plotL6Sweep(l6FigureRoot, weightTable, spectra, ...
    [-3 -1 0 0.2 0.4 0.6 0.8 1 1.3 2 5 10], ...
    fullfile(outputRoot, '3.1_L6_Weight_Sweep_Eigenvalue_Histograms.pdf'), ...
    binWidth);

plotSweep(weightTable, spectra, 'gammaI_sweep', 'wI', 'gammaI', ...
    [0.85 0.39 0.18], ...
    'Inhibition derivative-weight sweep', ...
    fullfile(outputRoot, '3.2_Inhibition_Weight_Sweep_Eigenvalue_Histograms.pdf'), ...
    binWidth);
end

function plotL6Sweep(figureRoot, weightTable, spectra, weights, outputFile, binWidth)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [30 30 1900 1400]);
layout = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
faceColor = [0.18 0.47 0.72];

for panelIndex = 1:numel(weights)
    weight = weights(panelIndex);
    if weight < 0
        canonicalIndex = find(strcmp(weightTable.sweepType, 'gamma6_sweep') & ...
            abs(weightTable.w6 - weight) < 1e-12, 1);
        if isempty(canonicalIndex)
            error('Missing canonical negative L6 spectrum for w=%g.', weight);
        end
        realPart = real(spectra{canonicalIndex}(:));
    else
        sourceFigure = fullfile(figureRoot, sprintf( ...
            'eigenvalue_realpart_hist_5Dmlph96baseline_L6w%.1f_angle_0.00_contr100.fig', ...
            weight));
        if ~isfile(sourceFigure)
            error('Missing L6 histogram source figure: %s', sourceFigure);
        end

        sourceHandle = openfig(sourceFigure, 'invisible');
        sourceHistogram = findobj(sourceHandle, 'Type', 'Histogram');
        if isempty(sourceHistogram)
            close(sourceHandle);
            error('No Histogram object found in %s.', sourceFigure);
        end
        realPart = sourceHistogram(1).Data(:);
        close(sourceHandle);
    end
    realPart = realPart(isfinite(realPart));
    if numel(realPart) ~= 4800
        error('Panel w=%g contains %d finite eigenvalues, expected 4800.', ...
            weight, numel(realPart));
    end

    edges = fixedWidthAdaptiveEdges(realPart, binWidth);
    ax = nexttile(layout);
    histogram(ax, realPart, edges, 'FaceColor', faceColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.92);
    hold(ax, 'on');
    xline(ax, 0, ':', 'Color', [0.20 0.20 0.20], 'LineWidth', 0.9);
    if edges(1) <= 1 && edges(end) >= 1
        xline(ax, 1, '--', 'Color', [0.75 0.12 0.12], 'LineWidth', 1.1);
    end
    xlim(ax, [edges(1) edges(end)]);
    grid(ax, 'on');
    ax.GridAlpha = 0.16;
    ax.FontSize = 10;
    xlabel(ax, 'Re(\lambda)');
    ylabel(ax, 'Count');
    title(ax, sprintf(['w=%g, L6 gain=1-w=%+.1f\n' ...
        'range=[%.3f, %.3f], max=%.3f'], ...
        weight, 1 - weight, min(realPart), max(realPart), max(realPart)), ...
        'FontSize', 10.5, 'FontWeight', 'normal');
end

title(layout, sprintf(['L6 equilibrium-weight sweep: full Jacobian, ' ...
    'orientation 0^{\\circ}, contrast 100 ' ...
    '(bin width = %.2f)'], binWidth), ...
    'FontSize', 15, 'FontWeight', 'bold');

exportgraphics(fig, outputFile, 'ContentType', 'vector');
close(fig);
fprintf('Saved %s\n', outputFile);
end

function plotSweep(weightTable, spectra, sweepType, weightName, gainName, ...
        faceColor, figureTitle, outputFile, binWidth)
indices = find(strcmp(weightTable.sweepType, sweepType));
weights = weightTable.(weightName)(indices);
gains = weightTable.(gainName)(indices);
[weights, order] = sort(weights);
indices = indices(order);
gains = gains(order);

if numel(indices) ~= 8
    error('Expected eight %s spectra, found %d.', sweepType, numel(indices));
end

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [50 50 1900 980]);
layout = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for panelIndex = 1:numel(indices)
    eigenvalues = spectra{indices(panelIndex)}(:);
    realPart = real(eigenvalues);
    realPart = realPart(isfinite(realPart));
    if numel(realPart) ~= 4800
        error('Panel w=%g contains %d finite eigenvalues, expected 4800.', ...
            weights(panelIndex), numel(realPart));
    end

    edges = fixedWidthAdaptiveEdges(realPart, binWidth);
    ax = nexttile(layout);
    histogram(ax, realPart, edges, 'FaceColor', faceColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.92);
    hold(ax, 'on');
    xline(ax, 0, ':', 'Color', [0.20 0.20 0.20], 'LineWidth', 0.9);
    if edges(1) <= 1 && edges(end) >= 1
        xline(ax, 1, '--', 'Color', [0.75 0.12 0.12], 'LineWidth', 1.1);
    end
    xlim(ax, [edges(1) edges(end)]);
    grid(ax, 'on');
    ax.GridAlpha = 0.16;
    ax.FontSize = 10;
    xlabel(ax, 'Re(\lambda)');
    ylabel(ax, 'Count');
    title(ax, sprintf(['w=%+.2f, gain=%+.2f\n' ...
        'range=[%.3f, %.3f], max=%.3f'], ...
        weights(panelIndex), gains(panelIndex), ...
        min(realPart), max(realPart), max(realPart)), ...
        'FontSize', 11, 'FontWeight', 'normal');
end

title(layout, sprintf(['%s: full Jacobian, orientation 0^{\\circ}, contrast 100 ' ...
    '(bin width = %.2f)'], ...
    figureTitle, binWidth), 'FontSize', 15, 'FontWeight', 'bold');

exportgraphics(fig, outputFile, 'ContentType', 'vector');
close(fig);
fprintf('Saved %s\n', outputFile);
end

function edges = fixedWidthAdaptiveEdges(values, binWidth)
low = floor(min(values) / binWidth) * binWidth - binWidth;
high = ceil(max(values) / binWidth) * binWidth + binWidth;
if high <= low
    high = low + binWidth;
end
edgeCount = round((high - low) / binWidth);
edges = low + (0:edgeCount) * binWidth;
if edges(end) < high - 10 * eps(max(abs(high), 1))
    edges(end + 1) = edges(end) + binWidth;
end
end
