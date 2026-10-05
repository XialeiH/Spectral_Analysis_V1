function result = run_ji_top_real_eigenspace_angles_task(taskIndex, sourceRoot, outputRoot)
% Plot the top-real clustered J_I eigenspace envelope for one angle.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('JI_SOURCE_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('JI_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrastValue = 100;
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > numel(angles)
    error('JIEigenspace:TaskIndex', 'Task index must be an integer from 1 to 4.');
end
if ~isfolder(sourceRoot)
    error('JIEigenspace:SourceRoot', 'Missing source root: %s', sourceRoot);
end
if isempty(outputRoot)
    error('JIEigenspace:OutputRoot', 'JI_OUTPUT_ROOT is required.');
end
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

angleValue = angles(taskIndex);
sourceFile = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    ['geometry_sections4_5_h96baseline_' ...
    'L6eqWS0p00_C0p00_I0p00_contr%d_angle_%.2f.mat'], ...
    contrastValue, angleValue));

result = plot_ji_top_real_eigenspace(sourceFile, outputRoot, 1e-3, ...
    angleValue, contrastValue);
end
