function plot_figure1g(timingFile, outputPdf)
% Plot paired timing distributions for the Paper 3 CG and h96 DNN surrogate.
timings = readtable(timingFile, 'FileType', 'text', 'Delimiter', '\t');
assert(height(timings) == 100, 'Figure 1G requires 100 timing trials.');

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172] / 255;
cgColor = neuroColors(1,:);
dnnColor = neuroColors(6,:);

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 16.5 7.3]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

draw_timing_panel(nexttile(layout), ...
    1000 .* timings.CGSecondsPerIteration, ...
    1000 .* timings.DNNSecondsPerIteration, ...
    'Average time per iteration (ms)', cgColor, dnnColor);
draw_timing_panel(nexttile(layout), ...
    timings.CGTotalSeconds, timings.DNNTotalSeconds, ...
    'End-to-end time for 50 iterations (s)', cgColor, dnnColor);

exportgraphics(fig, outputPdf, 'ContentType', 'vector');
end

function draw_timing_panel(ax, cgValues, dnnValues, yLabelText, cgColor, dnnColor)
hold(ax, 'on');
trialCount = numel(cgValues);
jitter = 0.16 .* sin((1:trialCount).' .* (sqrt(5)-1) .* pi);
scatter(ax, 1+jitter, cgValues, 38, cgColor, 'filled', ...
    'MarkerFaceAlpha', 0.32, 'MarkerEdgeAlpha', 0.32);
scatter(ax, 2+jitter, dnnValues, 38, dnnColor, 'filled', ...
    'MarkerFaceAlpha', 0.32, 'MarkerEdgeAlpha', 0.32);

[cgMedian, cgLow, cgHigh] = median_ci(cgValues);
[dnnMedian, dnnLow, dnnHigh] = median_ci(dnnValues);
errorbar(ax, 1, cgMedian, cgMedian-cgLow, cgHigh-cgMedian, ...
    'o', 'Color', cgColor, 'MarkerFaceColor', cgColor, ...
    'MarkerSize', 11, 'LineWidth', 3.0, 'CapSize', 18);
errorbar(ax, 2, dnnMedian, dnnMedian-dnnLow, dnnHigh-dnnMedian, ...
    'o', 'Color', dnnColor, 'MarkerFaceColor', dnnColor, ...
    'MarkerSize', 11, 'LineWidth', 3.0, 'CapSize', 18);

set(ax, 'YScale', 'log', 'XLim', [0.55 2.45], 'XTick', [1 2], ...
    'XTickLabel', {'CG', 'DNN surrogate'}, 'FontName', 'Arial', ...
    'FontSize', 26, 'LineWidth', 1.8, 'TickDir', 'out', 'Box', 'off');
ylabel(ax, yLabelText, 'FontSize', 29);
grid(ax, 'on');
ax.GridAlpha = 0.14;
speedup = median(cgValues ./ dnnValues);
text(ax, 1.5, 0.91, sprintf('Median speedup: %.1f\\times', speedup), ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontSize', 26, 'FontWeight', 'bold');
end

function [center, low, high] = median_ci(values)
values = sort(values(:));
center = median(values);
% Distribution-free 95% interval for the median at n=100.
low = values(40);
high = values(61);
end
