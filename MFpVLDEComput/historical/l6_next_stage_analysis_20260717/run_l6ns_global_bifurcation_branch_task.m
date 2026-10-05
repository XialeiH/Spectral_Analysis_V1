% Slurm-array entry point for the four signed secondary branches.

setupFile = getenv('L6NS_GLOBAL_BIFURCATION_SETUP');
if isempty(setupFile)
    error('L6NS_GLOBAL_BIFURCATION_SETUP is required.');
end
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_BRANCH_ROOT');
if isempty(outputDir)
    error('L6NS_GLOBAL_BIFURCATION_BRANCH_ROOT is required.');
end
taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId)
    taskId = str2double(getenv('L6NS_BRANCH_TASK_ID'));
end
if ~isfinite(taskId)
    error('A branch task ID in 1:4 is required.');
end
loaded = load(setupFile,'setup');
result = l6ns_global_bifurcation_branch(loaded.setup,taskId,outputDir); %#ok<NASGU>
