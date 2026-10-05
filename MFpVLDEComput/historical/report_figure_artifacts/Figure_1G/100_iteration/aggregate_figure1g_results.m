function aggregate_figure1g_results(runRoot)
files = dir(fullfile(runRoot, 'results', 'trial_*.tsv'));
assert(numel(files) == 100, 'Expected 100 trial files, found %d.', numel(files));
parts = cell(numel(files), 1);
for index = 1:numel(files)
    parts{index} = readtable(fullfile(files(index).folder, files(index).name), ...
        'FileType', 'text', 'Delimiter', '\t');
end
timings = sortrows(vertcat(parts{:}), 'Trial');
writetable(timings, fullfile(runRoot, 'figure1g_timings.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(runRoot, 'figure1g_timings.mat'), 'timings');
fprintf('Aggregated %d Figure 1G trials.\n', height(timings));
end
