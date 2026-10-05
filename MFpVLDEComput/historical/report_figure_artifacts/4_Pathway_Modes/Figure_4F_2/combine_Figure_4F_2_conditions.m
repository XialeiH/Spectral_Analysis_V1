function combine_Figure_4F_2_conditions(conditionRoot, outputFile)
% Combine all 16 conditions for Figure 4F.2.

artifactRoot = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(conditionRoot)
    conditionRoot = fullfile(artifactRoot, 'condition_data');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(artifactRoot, 'Figure_4F_2_visual_fit_data.mat');
end
files = dir(fullfile(conditionRoot, 'Figure_4F_2_condition_*.mat'));
if numel(files) ~= 16
    error('Figure4F2:Conditions', ...
        'Expected 16 condition files, found %d.', numel(files));
end

clusterFractions = zeros(16, 3, 2);
singularFractions = zeros(16, 3, 2);
clusterCoefficients = zeros(16, 3, 3);
singularCoefficients = zeros(16, 3, 3);
angles = zeros(16, 1);
contrasts = zeros(16, 1);
modeCounts = zeros(16, 3);
gramRcond = zeros(16, 1);
for index = 1:16
    data = load(fullfile(files(index).folder, files(index).name));
    row = data.taskIndex;
    clusterFractions(row, :, :) = data.clusterFractions;
    singularFractions(row, :, :) = data.singularFractions;
    clusterCoefficients(row, :, :) = data.clusterCoefficients;
    singularCoefficients(row, :, :) = data.singularCoefficients;
    angles(row) = data.angleDeg;
    contrasts(row) = data.contrast;
    modeCounts(row, :) = data.modeCounts;
    gramRcond(row) = data.gramRcond;
    jacobianNames = data.jacobianNames;
    metricDescription = data.metricDescription;
end

fractionNames = ["Visual manifold", "Off-manifold"];
coefficientNames = ["Log gain", "Orientation", "Log contrast"];
clusterFractionMean = squeeze(mean(clusterFractions, 1));
clusterFractionSD = squeeze(std(clusterFractions, 0, 1));
singularFractionMean = squeeze(mean(singularFractions, 1));
singularFractionSD = squeeze(std(singularFractions, 0, 1));
clusterCoefficientMean = squeeze(mean(clusterCoefficients, 1));
clusterCoefficientSD = squeeze(std(clusterCoefficients, 0, 1));
singularCoefficientMean = squeeze(mean(singularCoefficients, 1));
singularCoefficientSD = squeeze(std(singularCoefficients, 0, 1));

save(outputFile, 'angles', 'contrasts', 'jacobianNames', ...
    'fractionNames', 'coefficientNames', 'modeCounts', 'gramRcond', ...
    'metricDescription', 'clusterFractions', 'singularFractions', ...
    'clusterCoefficients', 'singularCoefficients', ...
    'clusterFractionMean', 'clusterFractionSD', ...
    'singularFractionMean', 'singularFractionSD', ...
    'clusterCoefficientMean', 'clusterCoefficientSD', ...
    'singularCoefficientMean', 'singularCoefficientSD', '-v7');

rows = cell(0, 9);
modeNames = ["Leading eigen-cluster", "Finite-time singular response"];
for modeIndex = 1:2
    if modeIndex == 1
        fractionData = clusterFractions;
        coefficientData = clusterCoefficients;
    else
        fractionData = singularFractions;
        coefficientData = singularCoefficients;
    end
    for jacobianIndex = 1:3
        for metricIndex = 1:5
            if metricIndex <= 2
                metricName = fractionNames(metricIndex);
                unit = "";
                values = squeeze(fractionData(:, jacobianIndex, metricIndex));
            else
                metricName = coefficientNames(metricIndex - 2);
                units = ["", "deg", ""];
                unit = units(metricIndex - 2);
                values = squeeze(coefficientData(:, jacobianIndex, metricIndex - 2));
            end
            rows(end+1, :) = {modeNames(modeIndex), ...
                jacobianNames(jacobianIndex), metricName, unit, ...
                mean(values), std(values), min(values), max(values), numel(values)}; %#ok<AGROW>
        end
    end
end
summaryTable = cell2table(rows, 'VariableNames', ...
    {'Mode', 'Jacobian', 'Metric', 'Unit', 'Mean', 'SD', ...
    'Minimum', 'Maximum', 'N'});
writetable(summaryTable, strrep(outputFile, '.mat', '.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
fprintf('Saved combined Figure 4F.2 data to %s.\n', outputFile);
end
