function extract_Figure_4H_geometry(sourceFile, outputFile)
% Extract the small geometry subset needed by Figure 4H.

source = load(sourceFile, 'Section4', 'GeometryMetrics');
featureNames = string(source.Section4.FeatureNames(:));
taskTangents = double(source.Section4.TaskTangents);
fixedPointMaps = source.GeometryMetrics.FixedPointMaps;

orientationIndex = find(featureNames == "orientation_deg", 1);
contrastIndex = find(featureNames == "log_contrast", 1);
assert(~isempty(orientationIndex) && ~isempty(contrastIndex), ...
    'Orientation and log-contrast task tangents are required.');

qOrientation = taskTangents(:, orientationIndex);
qLogContrast = taskTangents(:, contrastIndex);
save(outputFile, 'featureNames', 'qOrientation', 'qLogContrast', ...
    'fixedPointMaps', '-v7');
fprintf('Saved Figure 4H geometry subset to %s.\n', outputFile);
end
