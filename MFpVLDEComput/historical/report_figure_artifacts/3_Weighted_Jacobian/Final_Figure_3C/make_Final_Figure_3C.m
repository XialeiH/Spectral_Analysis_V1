function make_Final_Figure_3C()
% Plot the selected E-to-E kernels and spectrum histogram maps from Figure 3C.

projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project';
analysisRoot = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model', 'l6_next_stage_analysis_20260717');
controlRoot = fullfile(analysisRoot, ...
    'spatial_symmetry_zero_controls_20260811');
gkeRoot = fullfile(analysisRoot, 'gaussian_kernel_experiments_20260813');
resultRoot = fullfile(controlRoot, 'reequilibrated_all_controls_results');
setupFile = fullfile(analysisRoot, 'results_global_bifurcation_20260722', ...
    'global_bifurcation_setup.mat');
runtimeRoot = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model', 'geometry_analysis', ...
    'torch_h96_dedupe_array', 'Utils');
utilsRoot = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model', 'Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main', 'Utils');
outputPdf = fullfile(tempdir, 'Figure_3C_editable_source.pdf');
outputFig = fullfile(projectRoot, 'Spectral_Analysis_Report', ...
    'Figures', 'Final draft', 'Figure 3C.fig');

addpath(analysisRoot);
addpath(controlRoot);
addpath(gkeRoot);
addpath(utilsRoot);
addpath(runtimeRoot, '-begin');

setupData = load(setupFile, 'setup');
setup = setupData.setup;
extraSpectrumRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'Figure_3C', 'gke_cache', 'complete_spectra');

saved = {
    load(fullfile(resultRoot, 'baseline.mat'))
    load(fullfile(resultRoot, 'remove_L4_smoothing.mat'))
    load(fullfile(extraSpectrumRoot, 'shape_inscribed_triangle.mat'))
    load(fullfile(extraSpectrumRoot, 'shape_left_half.mat'))
    load(fullfile(extraSpectrumRoot, 'kernel_noise_cv0p30.mat'))
};
operators = cell(5, 1);
[operators{1}, ~] = l6ns_full_spatial_control_operators( ...
    setup.Context, 'baseline');
[operators{2}, ~] = l6ns_full_spatial_control_operators( ...
    setup.Context, 'remove_L4_smoothing');
[operators{3}, ~] = gke_build_operators( ...
    setup.Context, 'shape_inscribed_triangle');
[operators{4}, ~] = gke_build_operators( ...
    setup.Context, 'shape_left_half');
[operators{5}, ~] = gke_build_operators( ...
    setup.Context, 'kernel_noise_cv0p30');

labels = {'Baseline', 'Flattened', 'Triangle', 'Left Half', ...
    'Kernel Noise CV 0.30'};
for index = 1:5
    systems(index).Label = labels{index}; %#ok<AGROW>
    systems(index).KernelEE = local_centered_kernel( ...
        operators{index}.C_SS, [40 40]); %#ok<AGROW>
    systems(index).Eigenvalues = saved{index}.eigenvalues(:); %#ok<AGROW>
    assert(numel(systems(index).Eigenvalues) == 4800, ...
        'Each spectrum must contain 4,800 eigenvalues.');
end

% Preserve the original Figure 3C scales by deriving them from all five systems.
kernelStack = cat(3, systems.KernelEE);
positiveKernel = kernelStack(kernelStack > 0);
kernelFloor = max(max(kernelStack, [], 'all') * 1e-6, min(positiveKernel));
for index = 1:5
    systems(index).LogKernelEE = log10(max( ...
        systems(index).KernelEE, kernelFloor));
end
kernelLimits = [log10(kernelFloor) log10(max(kernelStack, [], 'all'))];

allEigenvalues = vertcat(systems.Eigenvalues);
binWidth = 0.05;
xLimits = [floor(min(real(allEigenvalues)) / binWidth) * binWidth, ...
    max(1.08, ceil(max(real(allEigenvalues)) / binWidth) * binWidth)];
binEdges = xLimits(1):binWidth:xLimits(2);
if binEdges(end) < xLimits(2)
    binEdges(end + 1) = xLimits(2);
end
histogramMaximum = 0;
for index = 1:5
    systems(index).HistogramCounts = histcounts( ...
        real(systems(index).Eigenvalues), binEdges);
    assert(sum(systems(index).HistogramCounts) == 4800, ...
        'Histogram counts must sum to 4,800.');
    histogramMaximum = max(histogramMaximum, ...
        max(systems(index).HistogramCounts));
end

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172
] / 255;

selected = [1 2 3 5];
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.2 0.2 15 24], 'PaperPositionMode', 'auto');
rowY = [0.765 0.535 0.305 0.075];
panelHeight = 0.180;
panelWidth = 0.288;
kernelX = 0.150;
distributionX = 0.515;

for row = 1:numel(selected)
    systemIndex = selected(row);

    kernelAx = axes(fig, 'Position', ...
        [kernelX rowY(row) panelWidth panelHeight]);
    imagesc(kernelAx, systems(systemIndex).LogKernelEE);
    axis(kernelAx, 'image');
    set(kernelAx, 'XTick', [], 'YTick', [], 'FontSize', 26, ...
        'LineWidth', 1.2, 'Box', 'on', 'CLim', kernelLimits);
    disableDefaultInteractivity(kernelAx);
    kernelAx.Toolbar.Visible = 'off';
    colormap(kernelAx, jet(256));
    if row == 1
        title(kernelAx, 'E-to-E Kernel', 'FontSize', 30, ...
            'FontWeight', 'normal');
        cb = colorbar(kernelAx, 'Position', ...
            [kernelX + panelWidth + 0.012, rowY(row), 0.018, panelHeight]);
        cb.FontSize = 26;
        cb.LineWidth = 1.0;
        cb.Title.String = 'log_{10}';
        cb.Title.Interpreter = 'tex';
        cb.Title.FontSize = 26;
    end

    labelAx = axes(fig, 'Position', [0.060 rowY(row) 0.055 panelHeight], ...
        'Visible', 'off', 'HitTest', 'off', 'PickableParts', 'none');
    text(labelAx, 0.5, 0.5, systems(systemIndex).Label, ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 27, ...
        'FontWeight', 'bold', 'Rotation', 90, 'Clipping', 'off');

    distributionAx = axes(fig, 'Position', ...
        [distributionX rowY(row) panelWidth panelHeight]);
    local_realpart_histogram(distributionAx, binEdges, ...
        systems(systemIndex).HistogramCounts, xLimits, ...
        histogramMaximum, neuroColors);
    if row < numel(selected)
        xlabel(distributionAx, '');
    end
    disableDefaultInteractivity(distributionAx);
    distributionAx.Toolbar.Visible = 'off';
    if row == 1
        title(distributionAx, {'Realpart of Eigenvalue', 'Distribution'}, ...
            'FontSize', 28, 'FontWeight', 'normal');
    end
end

savefig(fig, outputFig);
exportgraphics(fig, outputPdf, 'ContentType', 'image', 'Resolution', 600);
close(fig);
fprintf('Saved %s\n', outputPdf);
end

function kernel = local_centered_kernel(matrix, mapSize)
kernel = reshape(full(matrix(1, :)), mapSize);
kernel = circshift(kernel, [floor(mapSize(1) / 2) floor(mapSize(2) / 2)]);
kernel = abs(kernel);
kernel = kernel / max(sum(kernel(:)), eps);
end

function local_realpart_histogram(ax, binEdges, counts, xLimits, ...
        histogramMaximum, neuroColors)
binCenters = binEdges(1:end-1) + diff(binEdges) / 2;
hold(ax, 'on');
bar(ax, binCenters, counts, 1, 'BaseValue', 0, ...
    'FaceColor', [0.62 0.64 0.67], 'FaceAlpha', 0.50, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
xline(ax, 0, ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.8, ...
    'HandleVisibility', 'off');
xline(ax, 1, '--', 'Color', neuroColors(1, :), 'LineWidth', 2.2, ...
    'HandleVisibility', 'off');
xlabel(ax, 'Real part of eigenvalue', 'FontSize', 26);
xlim(ax, xLimits);
ylim(ax, [0 1.10 * histogramMaximum]);
set(ax, 'FontSize', 26, 'LineWidth', 1.2, 'Box', 'on', ...
    'Layer', 'top', 'Color', 'w', 'XGrid', 'off', 'YGrid', 'off', ...
    'YAxisLocation', 'right');
ylabel(ax, 'Count', 'FontSize', 26);
end
