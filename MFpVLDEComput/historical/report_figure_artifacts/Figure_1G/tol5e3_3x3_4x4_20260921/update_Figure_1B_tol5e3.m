function update_Figure_1B_tol5e3(resultDir, outputDir)
% Replace ONLY the 3x3 and 4x4 CG/DNN timings in manuscript Figure 1B.
% Trial Seconds includes iterations/convergence checks, not preparation/loading.
% The saved SNN timings and all larger-field points remain unchanged.
if nargin < 1
    resultDir = fullfile(fileparts(mfilename('fullpath')), 'results');
end
if nargin < 2
    outputDir = '/Users/xialeihuang/Desktop/Neuroscience_Project/Spectral_Analysis_Report/Figures/Draft';
end
files = [dir(fullfile(resultDir,'cg_dense_*HC_repeat_*.tsv')); ...
    dir(fullfile(resultDir,'dnn_fast_*HC_repeat_*.tsv'))];
assert(numel(files) == 400, 'Require all 400 new trial results before plotting.');
parts = cell(numel(files),1);
for k = 1:numel(files)
    parts{k} = readtable(fullfile(files(k).folder,files(k).name), ...
        'FileType','text','Delimiter','\t','TextType','string');
end
T = sortrows(vertcat(parts{:}),{'Method','FieldHC','Repeat'});
assert(all(T.Converged) && all(T.FinalResidual < 5e-3));
assert(all(T.Iterations <= 200) && all(T.Seconds > 0));
names = ["CG" "DNN surrogate"];
palette = [178 24 43; 33 102 172]/255;
means = zeros(2,2); sems = zeros(2,2);
for method = 1:2
    for field = 1:2
        rows = T.Method == names(method) & T.FieldHC == field+2;
        assert(sum(rows) == 100);
        assert(isequal(sort(T.Repeat(rows)),(1:100)'));
        means(field,method) = mean(T.Seconds(rows));
        sems(field,method) = std(T.Seconds(rows))/sqrt(sum(rows));
    end
end

fig = openfig(fullfile(outputDir,'Figure 1F.fig'),'invisible');
cleanup = onCleanup(@() close(fig));
fig.Units = 'inches';
fig.Position = [0.4 0.4 (51/36)*[11.2 9.1]];
ax = findall(fig,'Type','axes');
assert(numel(ax) == 1);
ax.YLabel.String = 'Wall-clock time (s)';
points = findall(ax,'Type','scatter');
originalY = arrayfun(@(p)p.YData,points,'UniformOutput',false);
snn = [];
updatedScatters = 0;
for k = 1:numel(points)
    x = points(k).XData;
    color = points(k).CData;
    assert(size(color,1) == 1 && all(x == x(1)));
    [distance,method] = min(sum((palette-color).^2,2));
    if distance < 1e-12 && ismember(x(1),[1 2])
        rows = T.Method == names(method) & T.FieldHC == x(1)+2;
        points(k).YData = T.Seconds(rows)';
        updatedScatters = updatedScatters+1;
    elseif sum((color-[254 227 145]/255).^2) < 1e-12
        snn = points(k).YData(:);
    end
end
assert(updatedScatters == 4 && numel(snn) == 5);
for k = 1:numel(points)
    distance = min(sum((palette-points(k).CData).^2,2));
    if points(k).XData(1) > 2 || distance >= 1e-12
        assert(isequal(points(k).YData,originalY{k}), 'Unrequested scatter changed.');
    end
end
lines = findall(ax,'Type','line');
updatedLines = 0;
for k = 1:numel(lines)
    [distance,method] = min(sum((palette-lines(k).Color).^2,2));
    if distance < 1e-12 && isequal(lines(k).XData,1:8)
        y = lines(k).YData;
        y(1:2) = means(:,method);
        lines(k).YData = y;
        updatedLines = updatedLines+1;
    end
end
assert(updatedLines == 2);
bars = findall(ax,'Type','errorbar');
updatedBars = 0;
for k = 1:numel(bars)
    [distance,method] = min(sum((palette-bars(k).Color).^2,2));
    if distance < 1e-12 && numel(bars(k).XData) == 8
        y = bars(k).YData; y(1:2) = means(:,method);
        low = bars(k).YNegativeDelta; low(1:2) = sems(:,method);
        high = bars(k).YPositiveDelta; high(1:2) = sems(:,method);
        bars(k).YData = y;
        bars(k).YNegativeDelta = low;
        bars(k).YPositiveDelta = high;
        updatedBars = updatedBars+1;
    end
end
assert(updatedBars == 2);
% Keep the existing range unless a new trial falls below its lower bound.
allSeconds = [];
allFields = []; allMethods = strings(0,1); allSamples = []; allTolerances = [];
fieldSizes = [3 4 6 8 10 20 30 40];
for k = 1:numel(points)
    values = points(k).YData(:); n = numel(values);
    field = fieldSizes(points(k).XData(1));
    [distance,method] = min(sum((palette-points(k).CData).^2,2));
    if distance < 1e-12
        name = names(method); tolerance = 1e-6;
        if ismember(field,[3 4]), tolerance = 5e-3; end
    else
        name = "SNN"; tolerance = NaN;
    end
    allSeconds = [allSeconds; values]; %#ok<AGROW>
    allFields = [allFields; repmat(field,n,1)]; %#ok<AGROW>
    allMethods = [allMethods; repmat(name,n,1)]; %#ok<AGROW>
    allSamples = [allSamples; (1:n)']; %#ok<AGROW>
    allTolerances = [allTolerances; repmat(tolerance,n,1)]; %#ok<AGROW>
end
plotted = table(allFields,allMethods,allSamples,allSeconds,allTolerances, ...
    'VariableNames',{'FieldHC','Method','SampleWithinGroup','Seconds','RelativeStepTolerance'});
plotted = sortrows(plotted,{'Method','FieldHC','SampleWithinGroup'});
assert(height(plotted) == 615);
writetable(plotted,fullfile(outputDir,'Figure 1B_only_iterations_data.csv'));
lowExponent = floor(log10(min(allSeconds)));
ax.YLim(1) = 10^lowExponent;
exponents = lowExponent:floor(log10(max(allSeconds)));
ax.YTick = 10.^exponents;
ax.YTickLabel = arrayfun(@(x)sprintf('$10^{%d}$',x),exponents,'UniformOutput',false);
drawnow;
savefig(fig,fullfile(outputDir,'Figure 1B_only_iterations.fig'));
exportgraphics(fig,fullfile(outputDir,'Figure 1B_only_iterations.pdf'), ...
    'ContentType','vector','BackgroundColor','white','Padding',16);
writetable(T,fullfile(outputDir,'Figure 1B_only_iterations_tol5e-3_trials.csv'));
for field = 1:2
    fprintf('Field %dx%d: CG mean %.9g s; DNN mean %.9g s; CG/DNN %.9gx\n', ...
        field+2,field+2,means(field,1),means(field,2),means(field,1)/means(field,2));
    for method = 1:2
        rows = T.FieldHC == field+2 & T.Method == names(method);
        fprintf('  %s: iterations %d-%d; mean iterations %.3f; residual %.6g-%.6g\n', ...
            names(method),min(T.Iterations(rows)),max(T.Iterations(rows)), ...
            mean(T.Iterations(rows)),min(T.FinalResidual(rows)),max(T.FinalResidual(rows)));
    end
end
fprintf('Unchanged 3x3 SNN mean %.9g s; SNN/CG (3x3) %.9gx\n', ...
    mean(snn),mean(snn)/means(1,1));
save(fullfile(outputDir,'Figure 1B_only_iterations_tol5e-3_summary.mat'), ...
    'T','means','sems','snn');
end
