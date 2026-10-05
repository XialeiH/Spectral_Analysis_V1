% Slurm-array entry point for the eight fixed-weight root-census probes.

setupFile = getenv('L6NS_GLOBAL_BIFURCATION_SETUP');
branchRoot = getenv('L6NS_GLOBAL_BIFURCATION_BRANCH_ROOT');
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_MULTISTART_ROOT');
if isempty(setupFile) || isempty(branchRoot) || isempty(outputDir)
    error('Global bifurcation setup, branch, and multistart roots are required.');
end
taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId)
    taskId = str2double(getenv('L6NS_MULTISTART_TASK_ID'));
end
if ~isfinite(taskId)
    error('A multistart task ID in 1:8 is required.');
end
loaded = load(setupFile,'setup');
result = l6ns_global_bifurcation_multistart( ...
    loaded.setup,branchRoot,taskId,outputDir); %#ok<NASGU>
