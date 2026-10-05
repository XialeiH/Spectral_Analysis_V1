function aggregate_figure1g_densemul(runRoot)
files = dir(fullfile(runRoot,'results_densemul','field_*_repeat_*.tsv'));
assert(~isempty(files),'No dense-multiplication result files found.');
parts = cell(numel(files),1);
for k = 1:numel(files)
    parts{k} = readtable(fullfile(files(k).folder,files(k).name), ...
        'FileType','text','Delimiter','\t');
end
timings = sortrows(vertcat(parts{:}),{'FieldHC','Repeat'});
trials4 = timings(timings.FieldHC==4,:);
assert(height(trials4)==100,'Dense 4x4 benchmark requires 100 trials.');
scaling = timings(timings.Repeat>=1 & timings.Repeat<=5,:);
assert(isequal(unique(scaling.FieldHC).',[4 6 8 10 20 30 40]), ...
    'Dense scaling benchmark is missing a field size.');
writetable(trials4,fullfile(runRoot,'figure1g_densemul_4x4_100trials.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(scaling,fullfile(runRoot,'field_scaling_densemul_timings.tsv'), ...
    'FileType','text','Delimiter','\t');
end
