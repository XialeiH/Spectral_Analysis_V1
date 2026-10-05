function combine_Figure_4F_3_conditions(conditionRoot, outputFile)
% Combine the 16 order-independent Shapley partitions for Figure 4F.3.

artifactRoot = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(conditionRoot)
    conditionRoot = fullfile(artifactRoot, 'condition_data');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(artifactRoot, 'Figure_4F_3_shapley_data.mat');
end
files = dir(fullfile(conditionRoot, 'Figure_4F_2_condition_*.mat'));
if numel(files) ~= 16
    error('Figure4F3:Conditions', ...
        'Expected 16 condition files, found %d.', numel(files));
end

clusterByCondition = zeros(16, 3, 4);
singularByCondition = zeros(16, 3, 4);
angles = zeros(16, 1);
contrasts = zeros(16, 1);
for index = 1:16
    data = load(fullfile(files(index).folder, files(index).name));
    row = data.taskIndex;
    clusterByCondition(row, :, :) = ...
        [data.clusterShapley, data.clusterFractions(:, 2)];
    singularByCondition(row, :, :) = ...
        [data.singularShapley, data.singularFractions(:, 2)];
    angles(row) = data.angleDeg;
    contrasts(row) = data.contrast;
    jacobianNames = data.jacobianNames;
    metricDescription = data.metricDescription;
end
if max(abs(sum(clusterByCondition, 3) - 1), [], 'all') > 1e-10 || ...
        max(abs(sum(singularByCondition, 3) - 1), [], 'all') > 1e-10
    error('Figure4F3:Partition', 'Shapley plus off-manifold does not sum to one.');
end

componentNames = ["Gain", "Orientation", "Contrast", "Off-manifold"];
clusterMean = squeeze(mean(clusterByCondition, 1));
clusterSD = squeeze(std(clusterByCondition, 0, 1));
singularMean = squeeze(mean(singularByCondition, 1));
singularSD = squeeze(std(singularByCondition, 0, 1));
save(outputFile, 'angles', 'contrasts', 'jacobianNames', ...
    'componentNames', 'metricDescription', 'clusterByCondition', ...
    'singularByCondition', 'clusterMean', 'clusterSD', ...
    'singularMean', 'singularSD', '-v7');

rows = cell(0, 8);
modeNames = ["Leading eigen-cluster", "Finite-time singular response"];
for modeIndex = 1:2
    if modeIndex == 1; values = clusterByCondition; else; values = singularByCondition; end
    for jacobian = 1:3
        for component = 1:4
            samples = squeeze(values(:, jacobian, component));
            rows(end+1, :) = {modeNames(modeIndex), jacobianNames(jacobian), ...
                componentNames(component), mean(samples), std(samples), ...
                min(samples), max(samples), numel(samples)}; %#ok<AGROW>
        end
    end
end
summaryTable = cell2table(rows, 'VariableNames', ...
    {'Mode', 'Jacobian', 'Component', 'Mean', 'SD', 'Minimum', 'Maximum', 'N'});
writetable(summaryTable, strrep(outputFile, '.mat', '.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
fprintf('Saved combined Figure 4F.3 data to %s.\n', outputFile);
end
