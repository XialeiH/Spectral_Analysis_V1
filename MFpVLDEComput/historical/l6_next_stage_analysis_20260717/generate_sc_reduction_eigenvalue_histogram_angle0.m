function generate_sc_reduction_eigenvalue_histogram_angle0(fullFile, reducedFile, outputRoot)
% Compare full and S/C-reduced spectra at angle 0 and contrast 100.

fullData = load(fullFile, 'eigenvalueData');
reducedData = load(reducedFile);
fullEigenvalues = fullData.eigenvalueData.J_baseline(:);
reducedEigenvalues = reducedData.reducedEigenvalues(:);

threshold = 0.05;
fullNearZero = nnz(abs(fullEigenvalues) <= threshold);
reducedNearZero = nnz(abs(reducedEigenvalues) <= threshold);
removedNearZero = fullNearZero - reducedNearZero;
removedPercent = 100 * removedNearZero / fullNearZero;
if fullNearZero ~= 2976 || reducedNearZero ~= 1369
    error('SCComparison:CountMismatch', ...
        'Unexpected near-zero counts: full=%d, reduced=%d.', ...
        fullNearZero, reducedNearZero);
end

binWidth = 0.05;
allReal = [real(fullEigenvalues); real(reducedEigenvalues)];
edgeMin = binWidth * floor(min(allReal) / binWidth);
edgeMax = binWidth * ceil(max(allReal) / binWidth);
edges = edgeMin:binWidth:edgeMax;
if numel(edges) < 2
    edges = edgeMin + [0 binWidth];
end
fullCounts = histcounts(real(fullEigenvalues), edges);
reducedCounts = histcounts(real(reducedEigenvalues), edges);
centers = edges(1:end-1) + binWidth / 2;

if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end
stem = fullfile(outputRoot, ...
    '10.2_Angle0_Full_vs_Reduced_Eigenvalue_Histogram');

fig = figure('Color', 'w', 'Position', [100 100 1160 700]);
set(fig, 'ToolBar', 'none');
hold on
stairs(edges, [fullCounts fullCounts(end)], 'LineWidth', 2.3, ...
    'Color', [0.0000 0.4470 0.7410], ...
    'DisplayName', 'Full S-C-I system (4,800 eigenvalues)');
stairs(edges, [reducedCounts reducedCounts(end)], 'LineWidth', 2.3, ...
    'Color', [0.8500 0.3250 0.0980], ...
    'DisplayName', 'Reduced E-I system (3,200 eigenvalues)');
xline(0, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2, ...
    'HandleVisibility', 'off');
hold off

xlim([edgeMin edgeMax]);
xlabel('Real part of eigenvalue, Re(\lambda)');
ylabel(sprintf('Eigenvalue count per %.2f-wide bin', binWidth));
title({'10.2  Eigenvalue histogram before and after S/C reduction', ...
    'Angle 0^\circ, contrast 100; common bin width 0.05'});
legend('Location', 'northwest');
grid on
box on
set(gca, 'FontSize', 14, 'LineWidth', 1.0);
disableDefaultInteractivity(gca);
axtoolbar(gca, {});

annotationText = sprintf([ ...
    'Exact near-zero disk, |\\lambda| \\leq %.2f\n' ...
    'Full: %s\nReduced: %s\nRemoved: %s (%.1f%%)'], ...
    threshold, local_integer(fullNearZero), local_integer(reducedNearZero), ...
    local_integer(removedNearZero), removedPercent);
annotation(fig, 'textbox', [0.15 0.19 0.27 0.22], ...
    'String', annotationText, 'Interpreter', 'tex', ...
    'FontSize', 13, 'BackgroundColor', 'white', ...
    'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 1.0, ...
    'FitBoxToText', 'on');

exportgraphics(fig, [stem '.pdf'], 'ContentType', 'vector');
savefig(fig, [stem '.fig']);

histogramTable = table(centers(:), fullCounts(:), reducedCounts(:), ...
    'VariableNames', {'realEigenvalueBinCenter', 'full4800Count', ...
    'reduced3200Count'});
writetable(histogramTable, [stem '_Data.tsv'], ...
    'FileType', 'text', 'Delimiter', '\t');

summaryTable = table(0, 100, numel(fullEigenvalues), ...
    numel(reducedEigenvalues), threshold, fullNearZero, reducedNearZero, ...
    removedNearZero, removedPercent, ...
    'VariableNames', {'angleDeg', 'contrast', 'fullDimension', ...
    'reducedDimension', 'nearZeroThreshold', 'fullNearZeroCount', ...
    'reducedNearZeroCount', 'removedNearZeroCount', ...
    'removedPercentOfFullNearZero'});
writetable(summaryTable, [stem '_Summary.tsv'], ...
    'FileType', 'text', 'Delimiter', '\t');
end

function value = local_integer(number)
value = sprintf('%.0f', number);
for index = (numel(value) - 2):-3:2
    value = [value(1:index-1) ',' value(index:end)]; %#ok<AGROW>
end
end
