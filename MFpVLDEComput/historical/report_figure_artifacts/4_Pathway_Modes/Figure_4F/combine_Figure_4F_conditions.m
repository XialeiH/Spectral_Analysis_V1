function combine_Figure_4F_conditions(conditionRoot, outputFile)
% Combine the 16 Figure 4F condition outputs.

if nargin < 1 || isempty(conditionRoot)
    conditionRoot = fullfile(fileparts(mfilename('fullpath')), 'condition_data');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(fileparts(mfilename('fullpath')), ...
        'Figure_4F_projection_data.mat');
end
files = dir(fullfile(conditionRoot, 'Figure_4F_condition_*.mat'));
if numel(files) ~= 16
    error('Figure4F:Conditions', 'Expected 16 condition files, found %d.', numel(files));
end

clusterByCondition = zeros(16, 3, 4);
singularByCondition = zeros(16, 3, 4);
angles = zeros(16, 1);
contrasts = zeros(16, 1);
modeCounts = zeros(16, 3);
for index = 1:16
    data = load(fullfile(files(index).folder, files(index).name));
    row = data.taskIndex;
    clusterByCondition(row, :, :) = data.clusterProjections;
    singularByCondition(row, :, :) = data.singularProjections;
    angles(row) = data.angleDeg;
    contrasts(row) = data.contrast;
    modeCounts(row, :) = data.modeCounts;
    jacobianNames = data.jacobianNames;
end
componentNames = ["Gain", "Orientation", "Contrast", "Orthogonal residual"];
clusterMean = squeeze(mean(clusterByCondition, 1));
clusterSD = squeeze(std(clusterByCondition, 0, 1));
singularMean = squeeze(mean(singularByCondition, 1));
singularSD = squeeze(std(singularByCondition, 0, 1));
save(outputFile, 'angles', 'contrasts', 'jacobianNames', 'componentNames', ...
    'modeCounts', 'clusterByCondition', 'singularByCondition', ...
    'clusterMean', 'clusterSD', 'singularMean', 'singularSD', '-v7');

rows = cell(0, 8);
for modeIndex = 1:2
    if modeIndex == 1
        modeName = "Leading eigen-cluster";
        means = clusterMean;
        deviations = clusterSD;
    else
        modeName = "Finite-time singular response";
        means = singularMean;
        deviations = singularSD;
    end
    for jacobianIndex = 1:3
        for componentIndex = 1:4
            rows(end+1, :) = {modeName, jacobianNames(jacobianIndex), ...
                componentNames(componentIndex), means(jacobianIndex, componentIndex), ...
                deviations(jacobianIndex, componentIndex), 16, ...
                min(squeeze((modeIndex == 1) * clusterByCondition(:, jacobianIndex, componentIndex) + ...
                    (modeIndex == 2) * singularByCondition(:, jacobianIndex, componentIndex))), ...
                max(squeeze((modeIndex == 1) * clusterByCondition(:, jacobianIndex, componentIndex) + ...
                    (modeIndex == 2) * singularByCondition(:, jacobianIndex, componentIndex)))}; %#ok<AGROW>
        end
    end
end
summaryTable = cell2table(rows, 'VariableNames', ...
    {'Mode', 'Jacobian', 'Component', 'Mean', 'SD', 'N', 'Minimum', 'Maximum'});
writetable(summaryTable, strrep(outputFile, '.mat', '.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
fprintf('Saved combined Figure 4F data to %s.\n', outputFile);
end
