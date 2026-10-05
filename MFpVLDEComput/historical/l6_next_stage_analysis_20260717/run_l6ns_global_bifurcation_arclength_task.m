% Slurm-array entry point for pseudo-arclength continuation of four branches.

setupFile = getenv('L6NS_GLOBAL_BIFURCATION_SETUP');
refinedRoot = getenv('L6NS_GLOBAL_BIFURCATION_REFINED_ROOT');
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_ARCLENGTH_ROOT');
if isempty(setupFile) || isempty(refinedRoot) || isempty(outputDir)
    error('Setup, refined branch, and arclength branch roots are required.');
end
taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId)
    taskId = str2double(getenv('L6NS_BRANCH_TASK_ID'));
end
if ~isfinite(taskId) || taskId<1 || taskId>4
    error('An arclength task ID in 1:4 is required.');
end
names = {'l6','l6','inhibition','inhibition'};
signs = [-1 1 -1 1];
stem = sprintf('%s_sign_%+d_result.mat',names{taskId},signs(taskId));
stem = strrep(stem,'+','p');
stem = strrep(stem,'-','m');
setupLoaded = load(setupFile,'setup');
branchLoaded = load(fullfile(refinedRoot,stem),'result');
result = l6ns_global_bifurcation_arclength(setupLoaded.setup, ...
    branchLoaded.result,outputDir); %#ok<NASGU>
