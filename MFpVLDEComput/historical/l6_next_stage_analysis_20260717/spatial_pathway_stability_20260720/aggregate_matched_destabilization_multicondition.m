function summary = aggregate_matched_destabilization_multicondition(outputRoot)
% Aggregate the 16 matched-destabilization conditions and plot diagnostics.

if nargin < 1 || isempty(outputRoot)
    outputRoot = getenv('MATCHED_OUTPUT_ROOT');
end
dataRoot = fullfile(outputRoot, 'data');
figureRoot = fullfile(outputRoot, 'figures');
files = dir(fullfile(dataRoot, 'matched_destabilization_Angle_*_Contrast_*.tsv'));
assert(numel(files) == 16, 'Expected 16 condition summaries, found %d.', numel(files));

tables = cell(numel(files), 1);
for k = 1:numel(files)
    tables{k} = readtable(fullfile(files(k).folder, files(k).name), ...
        'FileType', 'text', 'Delimiter', '\t');
end
summary = vertcat(tables{:});
summary = sortrows(summary, {'contrast', 'angle'});
writetable(summary, fullfile(dataRoot, 'matched_destabilization_all_conditions.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(dataRoot, 'matched_destabilization_all_conditions.mat'), ...
    'summary', '-v7.3');

angles = unique(summary.angle).';
contrasts = unique(summary.contrast).';
metrics = {'rightModeOverlap','leftModeOverlap','trackedGapI', ...
    'alphaMatchError','gammaI','cStab'};
titles = {'Tracked right-mode overlap', 'Tracked left-mode overlap', ...
    'I-path tracked-to-top branch gap', 'Full-J spectral-abscissa match error', ...
    'Matched inhibitory gain  \gamma_I', 'First-order compensation  c_{stab}'};

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [30 30 1550 1080]);
t = tiledlayout(fig, 2, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
for metricIndex = 1:numel(metrics)
    values = nan(numel(contrasts), numel(angles));
    for contrastIndex = 1:numel(contrasts)
        for angleIndex = 1:numel(angles)
            row = summary.contrast == contrasts(contrastIndex) & ...
                abs(summary.angle - angles(angleIndex)) < 1e-9;
            values(contrastIndex, angleIndex) = summary.(metrics{metricIndex})(row);
        end
    end
    ax = nexttile(t);
    imagesc(ax, angles, contrasts, values);
    set(ax, 'YDir', 'normal');
    colormap(ax, parula(256));
    colorbar(ax);
    if metricIndex > 3
        xlabel(ax, 'orientation (deg)');
    end
    ylabel(ax, 'contrast');
    title(ax, titles{metricIndex}, 'FontWeight', 'bold');
    xticks(ax, angles); yticks(ax, contrasts);
    for contrastIndex = 1:numel(contrasts)
        for angleIndex = 1:numel(angles)
            text(ax, angles(angleIndex), contrasts(contrastIndex), ...
                sprintf('%.3g', values(contrastIndex, angleIndex)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, ...
                'Color', contrast_text_color(values, values(contrastIndex, angleIndex)));
        end
    end
end
title(t, {['Matched L6-increase versus inhibition-decrease mechanism across ', ...
    '16 audited conditions'], ...
    ['Each pair is matched by the full-J spectral abscissa; mode overlaps use ', ...
    'the same baseline-continuation branch']}, ...
    'FontSize', 15, 'FontWeight', 'bold');
set(fig, 'Renderer', 'painters');
summaryPdf = fullfile(figureRoot, ...
    'Matched_Destabilization_Across_Conditions_Quantitative_Summary.pdf');
exportgraphics(fig, summaryPdf, 'ContentType', 'vector');
close(fig);

assert(all(summary.reconstructionError < 1e-12), ...
    'At least one pathway reconstruction check failed.');
fprintf(['Aggregated %d conditions. Median right overlap %.6f; median left ', ...
    'overlap %.6f; max alpha match error %.3e.\n'], height(summary), ...
    median(summary.rightModeOverlap), median(summary.leftModeOverlap), ...
    max(summary.alphaMatchError));
disp(summary(:, {'angle','contrast','alphaL6','alphaI','trackedGapI', ...
    'rightModeOverlap','leftModeOverlap','minimumRightSubspaceCosine', ...
    'minimumLeftSubspaceCosine'}));
end

function color = contrast_text_color(values, value)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    color = [0 0 0];
    return
end
fraction = (value - min(finiteValues)) / ...
    max(max(finiteValues) - min(finiteValues), eps);
if fraction < 0.55
    color = [1 1 1];
else
    color = [0 0 0];
end
end
