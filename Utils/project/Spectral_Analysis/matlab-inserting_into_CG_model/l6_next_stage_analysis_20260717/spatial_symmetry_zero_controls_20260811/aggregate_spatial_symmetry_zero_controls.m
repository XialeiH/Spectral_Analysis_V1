function summary = aggregate_spatial_symmetry_zero_controls(resultRoots, outputRoot)
% Aggregate repeated spatial-permutation zero-count controls.

if nargin < 1 || isempty(resultRoots)
    baseRoot = repro_paths('project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717/spatial_symmetry_zero_controls_20260811');
    resultRoots = [{fullfile(baseRoot, 'results')}, ...
        arrayfun(@(seed)fullfile(baseRoot, 'results_mc', sprintf('seed%d', seed)), ...
        8112027:8112030, 'UniformOutput', false)];
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = fullfile(repro_paths('project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717/spatial_symmetry_zero_controls_20260811'), 'results_mc');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

allRows = table();
for replicate = 1:numel(resultRoots)
    source = fullfile(resultRoots{replicate}, 'zero_count_summary.tsv');
    if ~isfile(source)
        error('Missing replicate summary: %s', source);
    end
    rows = readtable(source, 'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string');
    rows.replicate = repmat(replicate, height(rows), 1);
    allRows = [allRows; rows]; %#ok<AGROW>
end

controls = unique(allRows.control, 'stable');
summaryRows = cell(numel(controls), 16);
for index = 1:numel(controls)
    rows = allRows(allRows.control == controls(index), :);
    summaryRows(index, :) = {controls(index), height(rows), ...
        mean(rows.countAbsEigLE1e2), min(rows.countAbsEigLE1e2), ...
        max(rows.countAbsEigLE1e2), ...
        mean(rows.countAbsEigLE2p5e2), min(rows.countAbsEigLE2p5e2), ...
        max(rows.countAbsEigLE2p5e2), ...
        mean(rows.countAbsEigLE5e2), min(rows.countAbsEigLE5e2), ...
        max(rows.countAbsEigLE5e2), ...
        mean(rows.countAbsEigLE1e10), min(rows.countAbsEigLE1e10), ...
        max(rows.countAbsEigLE1e10), ...
        mean(rows.maximumRealEigenvalue), std(rows.maximumRealEigenvalue)};
end
summary = cell2table(summaryRows, 'VariableNames', ...
    {'control', 'replicateCount', 'meanN001', 'minimumN001', 'maximumN001', ...
    'meanN0025', 'minimumN0025', 'maximumN0025', ...
    'meanN005', 'minimumN005', 'maximumN005', ...
    'meanExactZero', 'minimumExactZero', 'maximumExactZero', ...
    'meanMaximumRealEigenvalue', 'stdMaximumRealEigenvalue'});

writetable(allRows, fullfile(outputRoot, 'all_replicate_zero_counts.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
writetable(summary, fullfile(outputRoot, 'monte_carlo_zero_count_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

selected = ~ismember(summary.control, ...
    ["baseline", "shuffle_values_preserve_zero_mask", "block_mean_only"]);
plotRows = summary(selected, :);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1450 800]);
errorbar(1:height(plotRows), plotRows.meanN005, ...
    plotRows.meanN005 - plotRows.minimumN005, ...
    plotRows.maximumN005 - plotRows.meanN005, ...
    'o', 'LineWidth', 1.5, 'MarkerFaceColor', [0.18 0.48 0.72]);
hold on;
yline(2976, 'k--', 'Baseline N_{0.05}=2976', 'LineWidth', 1.2);
grid on; box on;
xticks(1:height(plotRows));
xticklabels(local_display_labels(plotRows.control));
xtickangle(22);
ylabel('N_{0.05}');
title('Near-zero counts across five targeted spatial-permutation realizations');
exportgraphics(fig, fullfile(outputRoot, ...
    '04_monte_carlo_near_zero_count_ranges.pdf'), 'ContentType', 'vector');
close(fig);
end

function labels = local_display_labels(controls)
labels = strings(size(controls));
for index = 1:numel(controls)
    switch controls(index)
        case "preserve_translation_break_reflection_C2"
            labels(index) = "T kept; C2 broken";
        case "preserve_reflection_C2_break_translation"
            labels(index) = "C2 kept; T broken";
        case "break_translation_and_reflection_shared_SC"
            labels(index) = "T/C2 broken; S-C kept";
        case "break_both_and_SC_alignment"
            labels(index) = "T/C2/S-C broken";
        otherwise
            labels(index) = replace(controls(index), '_', ' ');
    end
end
end
