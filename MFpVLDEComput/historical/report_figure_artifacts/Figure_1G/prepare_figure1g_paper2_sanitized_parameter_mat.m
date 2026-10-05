function prepare_figure1g_paper2_sanitized_parameter_mat(sourceFile, outputFile)
% Remove workstation-only state while preserving the Paper 2 model inputs.

sourceInfo = whos('-file', sourceFile);
payload = load(sourceFile);

staleNames = { ...
    'CurrentFolder', 'DataFolder', 'DataFolder1', 'DriveDataFolder', ...
    'SaveFolder', 'FigurePath', 'FigurePaperPath', ...
    'hFrVAll', 'hFrVOnOff'};
removedNames = intersect(fieldnames(payload), staleNames, 'stable');
if ~isempty(removedNames)
    payload = rmfield(payload, removedNames);
end

outputDir = fileparts(outputFile);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
save(outputFile, '-struct', 'payload', '-v7');

outputInfo = whos('-file', outputFile);
expectedNames = setdiff({sourceInfo.name}, removedNames, 'stable');
actualNames = {outputInfo.name};
assert(isequal(sort(expectedNames), sort(actualNames)), ...
    'Sanitized MAT variable names do not match the retained source variables.');

for idx = 1:numel(expectedNames)
    name = expectedNames{idx};
    sourceVar = sourceInfo(strcmp({sourceInfo.name}, name));
    outputVar = outputInfo(strcmp(actualNames, name));
    assert(isequal(sourceVar.size, outputVar.size), ...
        'Size mismatch for retained variable %s.', name);
    assert(strcmp(sourceVar.class, outputVar.class), ...
        'Class mismatch for retained variable %s.', name);
end

fprintf('Sanitized %s -> %s\n', sourceFile, outputFile);
fprintf('Removed staging-only variables: %s\n', strjoin(removedNames, ', '));
fprintf('Verified %d retained variables by name, class, and size.\n', numel(expectedNames));
end
