function aggregate_figure1g_scaling(runRoot)
files = dir(fullfile(runRoot,'results','field_*_repeat_*.tsv'));
assert(~isempty(files),'No field-scaling result files found.');
parts = cell(numel(files),1);
for k = 1:numel(files)
    parts{k} = readtable(fullfile(files(k).folder,files(k).name), ...
        'FileType','text','Delimiter','\t');
end
timings = vertcat(parts{:});
timings = sortrows(timings,{'FieldHC','Repeat'});
writetable(timings,fullfile(runRoot,'field_scaling_timings.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(runRoot,'field_scaling_timings.mat'),'timings');
end
