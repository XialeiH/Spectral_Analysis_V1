function verify_published_numbers()
% Check retained data against regenerated graph objects, including hidden handles.
root=fileparts(fileparts(mfilename('fullpath')));
inputRoot=fullfile(root,'Data','project','Spectral_Analysis_Report','Figures','Draft');
expected=readtable(fullfile(inputRoot,'Figure 1B_only_iterations_data.csv'));
actual=readtable(fullfile(root,'Outputs','Main_01_B','Main_01_B_data.csv'));
assert(height(expected)==615 && height(actual)==615);
assert(isequal(expected.Seconds,actual.Seconds));
f=openfig(fullfile(root,'Outputs','Main_01_B','Main_01_B.fig'),'invisible');
scatters=findall(f,'Type','Scatter');
values=[];
for k=1:numel(scatters), values=[values;scatters(k).YData(:)]; end %#ok<AGROW>
assert(numel(values)==615);
assert(isequal(sort(values),sort(actual.Seconds)));
close(f);
f=openfig(fullfile(root,'Outputs','Main_01_C','Main_01_C.fig'),'invisible');
points=findall(f,'Type','Scatter');
assert(sum(arrayfun(@(p)numel(p.XData),points))==28000);
close(f);
reference=openfig(fullfile(inputRoot,'Figure 1E.fig'),'invisible');
generated=openfig(fullfile(root,'Outputs','Main_01_FG','Main_01_FG.fig'),'invisible');
before=findall(reference,'Type','image'); after=findall(generated,'Type','image');
assert(numel(before)==6 && numel(after)==6);
for k=1:6, assert(isequal(before(k).CData,after(k).CData)); end
close(reference); close(generated);
summary=readtable(fullfile(root,'Outputs','Main_04_B','Main_04_B_summary.csv'));
assert(isequal(summary.SmallSingularCount,[2981;1792]));
assert(isequal(summary.EffectiveRank,[1819;3008]));
assert(isequal(summary.NearZeroRealCount,[2988;1806]));
trajectoryRoot=fullfile('project','Spectral_Analysis','matlab-inserting_into_CG_model', ...
    'report_figure_artifacts','6_Pathway_Compensation','Figure_5D','fixed_input_dynamic_difference');
expectedPeaks=readtable(fullfile(root,'Data',trajectoryRoot,'fixed_input_peak_time.tsv'), ...
    'FileType','text','Delimiter','\t');
actualPeaks=readtable(fullfile(root,'Outputs','Main_06_D','fixed_input_peak_time.tsv'), ...
    'FileType','text','Delimiter','\t');
assert(height(expectedPeaks)==22 && height(actualPeaks)==22);
assert(isequal(expectedPeaks.fixedInputPeakTimeMs,actualPeaks.fixedInputPeakTimeMs));
assert(max(abs(expectedPeaks.fixedInputPeakHC-actualPeaks.fixedInputPeakHC))<1e-12);
fprintf('PASS: 615 timings, 28000 agreement points, six image arrays, L4 spectral counts, and 22 trajectory peaks.\n');
end
