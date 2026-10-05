function make_Final_Figure_2C
% Reformat Figure 2A.1 with one enlarged shared w=0 baseline panel.

artifactDir = fileparts(mfilename('fullpath'));
modelRoot = ['/Users/xialeihuang/Desktop/Neuroscience_Project/Spectral_Analysis/' ...
    'matlab-inserting_into_CG_model/Complete_Code_for_Paper3/' ...
    'NYU-Vision-2Drive-main'];
l6FigureDir = fullfile(modelRoot, 'Figures', ...
    'spectral_analysis_eigenvalue_eigenvectors', 'eigenspectrum_L6weight');
outputDir = ['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis_Report/Figures/Final draft'];
outputPdf = fullfile(outputDir, 'Figure 2C.pdf');
outputFig = fullfile(outputDir, 'Figure 2C.fig');

weights = [-1 0 1 5];
l6Spectra = cell(size(weights));
iSpectra = cell(size(weights));
for index = 1:numel(weights)
    weight = weights(index);
    if weight == -1
        l6Spectra{index} = load_cached_spectrum(artifactDir, 'L6', weight);
    else
        sourceFigure = fullfile(l6FigureDir, sprintf( ...
            'eigenvalue_spectrum_5Dmlph96baseline_L6w%.1f_angle_0.00_contr100.fig', ...
            weight));
        l6Spectra{index} = load_spectrum_from_figure(sourceFigure);
    end
    if weight == 0
        iSpectra{index} = l6Spectra{index};
    else
        iSpectra{index} = load_cached_spectrum(artifactDir, 'I', weight);
    end
end

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172] / 255;
l6Color = neuroColors(6,:);
iColor = neuroColors(1,:);
pointColor = [0.34 0.34 0.34];
axisColor = [0.12 0.12 0.12];

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.2 0.2 28 13], 'Renderer', 'opengl');

smallWidth = 0.145;
smallHeight = 0.365;
bottomY = 0.115;
topY = 0.555;
fullHeight = topY + smallHeight - bottomY;
xSmall = [0.065 0.615 0.805];
centerPosition = [0.245 bottomY 0.325 fullHeight];

% Left stacked pair: w=-1.
ax = axes(fig, 'Position', [xSmall(1) topY smallWidth smallHeight]);
plot_spectrum_histogram(ax, l6Spectra{1}, l6Color, pointColor, axisColor);
title(ax, '$w_6=-1$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel(ax, 'Imaginary part of eigenvalue', 'FontSize', 22);

ax = axes(fig, 'Position', [xSmall(1) bottomY smallWidth smallHeight]);
plot_spectrum_histogram(ax, iSpectra{1}, iColor, pointColor, axisColor);
title(ax, '$w_I=-1$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel(ax, 'Imaginary part of eigenvalue', 'FontSize', 22);

% The two original w=0 panels are identical; retain one at twice scale.
ax = axes(fig, 'Position', centerPosition);
plot_spectrum_histogram(ax, l6Spectra{2}, l6Color, pointColor, axisColor);
title(ax, '$w=0$ (Baseline)', 'Interpreter', 'latex', 'FontSize', 28);

% Right stacked pairs: w=1 and w=5.
rightWeights = [1 5];
for column = 1:2
    spectrumIndex = column + 2;
    ax = axes(fig, 'Position', ...
        [xSmall(column + 1) topY smallWidth smallHeight]);
    plot_spectrum_histogram(ax, l6Spectra{spectrumIndex}, ...
        l6Color, pointColor, axisColor);
    title(ax, sprintf('$w_6=%g$', rightWeights(column)), ...
        'Interpreter', 'latex', 'FontSize', 24);

    ax = axes(fig, 'Position', ...
        [xSmall(column + 1) bottomY smallWidth smallHeight]);
    maximumCount = plot_spectrum_histogram(ax, ...
        iSpectra{spectrumIndex}, iColor, pointColor, axisColor);
    title(ax, sprintf('$w_I=%g$', rightWeights(column)), ...
        'Interpreter', 'latex', 'FontSize', 24);
    if column == 2
        add_count_axis(fig, ax, maximumCount, axisColor);
    end
end

% Count scale for the upper-right panel.
upperRight = findobj(fig, 'Type', 'axes', 'Position', ...
    [xSmall(3) topY smallWidth smallHeight]);
maximumCount = max_histogram_count(l6Spectra{4});
add_count_axis(fig, upperRight(1), maximumCount, axisColor);

annotation(fig, 'textbox', [0.39 0.018 0.25 0.055], ...
    'String', 'Real part of eigenvalue', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontName', 'Arial', 'FontSize', 26, 'Color', axisColor);

savefig(fig, outputFig);
if isfile(outputPdf)
    delete(outputPdf);
end
exportgraphics(fig, outputPdf, 'ContentType', 'image', 'Resolution', 600, ...
    'BackgroundColor', 'white');
close(fig);
fprintf('Saved %s and %s.\n', outputPdf, outputFig);
end

function lambda = load_cached_spectrum(artifactDir, pathwayName, weight)
cacheFile = fullfile(artifactDir, sprintf('%s_weight_%+.1f_eigenvalues.mat', ...
    pathwayName, weight));
assert(isfile(cacheFile), 'Missing eigenspectrum cache: %s', cacheFile);
loaded = load(cacheFile, 'lambda');
lambda = loaded.lambda(:);
assert(numel(lambda) == 4800, 'Expected 4800 eigenvalues in %s.', cacheFile);
end

function lambda = load_spectrum_from_figure(figureFile)
assert(isfile(figureFile), 'Missing retained MATLAB figure: %s', figureFile);
sourceFigure = openfig(figureFile, 'invisible');
cleanup = onCleanup(@() close(sourceFigure));
objects = [findall(sourceFigure, 'Type', 'scatter'); ...
    findall(sourceFigure, 'Type', 'line')];
bestCount = 0;
bestX = [];
bestY = [];
for index = 1:numel(objects)
    try
        x = objects(index).XData;
        y = objects(index).YData;
    catch
        continue;
    end
    if isnumeric(x) && isnumeric(y) && numel(x) == numel(y) && ...
            numel(x) > bestCount
        bestCount = numel(x);
        bestX = x;
        bestY = y;
    end
end
assert(bestCount == 4800, ...
    'Could not identify the 4800-point eigenspectrum in %s.', figureFile);
lambda = bestX(:) + 1i * bestY(:);
end

function maximumCount = plot_spectrum_histogram(ax, lambda, ...
        histogramColor, pointColor, axisColor)
realPart = real(lambda(:));
imaginaryPart = imag(lambda(:));
imaginaryLimit = max(0.05, 1.08 * max(abs(imaginaryPart)));
minimumReal = min(realPart);
maximumReal = max(realPart);
leftLimit = minimumReal - 0.04 * max(1, abs(minimumReal));
rightLimit = max(1.10, maximumReal + 0.04 * max(1, abs(maximumReal)));

[counts, edges] = histcounts(realPart, 'BinWidth', 0.05);
centers = (edges(1:end-1) + edges(2:end)) / 2;
maximumCount = max(1, max(counts));
barTops = -imaginaryLimit + ...
    2 * imaginaryLimit * 0.88 * counts / maximumCount;

hold(ax, 'on');
bar(ax, centers, barTops, 1, 'BaseValue', -imaginaryLimit, ...
    'FaceColor', histogramColor, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.68, 'ShowBaseLine', 'off');
scatter(ax, realPart, imaginaryPart, 12, pointColor, 'filled', ...
    'MarkerFaceAlpha', 0.42, 'MarkerEdgeAlpha', 0.42);
xline(ax, 0, ':', 'Color', [0.38 0.38 0.38], 'LineWidth', 1.4);
xline(ax, 1, '--', 'Color', [227 74 51] / 255, 'LineWidth', 1.8);
xlim(ax, [leftLimit rightLimit]);
ylim(ax, [-imaginaryLimit imaginaryLimit]);
set(ax, 'FontName', 'Arial', 'FontSize', 20, 'LineWidth', 1.2, ...
    'Box', 'on', 'Layer', 'top', 'TickDir', 'out', 'Color', 'none', ...
    'XColor', axisColor, 'YColor', axisColor);
ax.Toolbar.Visible = 'off';
grid(ax, 'off');
end

function maximumCount = max_histogram_count(lambda)
counts = histcounts(real(lambda(:)), 'BinWidth', 0.05);
maximumCount = max(1, max(counts));
end

function add_count_axis(fig, sourceAxes, maximumCount, axisColor)
countAxes = axes(fig, 'Units', sourceAxes.Units, ...
    'Position', sourceAxes.Position, 'Color', 'none', ...
    'YAxisLocation', 'right', 'YLim', [0 maximumCount / 0.88], ...
    'XLim', sourceAxes.XLim, 'XTick', [], 'Box', 'off', ...
    'FontName', 'Arial', 'FontSize', 18, 'LineWidth', 1.2, ...
    'HitTest', 'off', 'YColor', axisColor);
countAxes.XAxis.Visible = 'off';
countAxes.Toolbar.Visible = 'off';
ylabel(countAxes, 'Count', 'FontName', 'Arial', 'FontSize', 22, ...
    'Color', axisColor);
end
