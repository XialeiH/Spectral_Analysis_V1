% Slurm-array entry point for adaptive continuation of four signed branches.

setupFile = getenv('L6NS_GLOBAL_BIFURCATION_SETUP');
coarseRoot = getenv('L6NS_GLOBAL_BIFURCATION_BRANCH_ROOT');
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_REFINED_ROOT');
if isempty(setupFile) || isempty(coarseRoot) || isempty(outputDir)
    error('Setup, coarse branch, and refined branch roots are required.');
end
taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId)
    taskId = str2double(getenv('L6NS_BRANCH_TASK_ID'));
end
if ~isfinite(taskId) || taskId<1 || taskId>4
    error('A refinement task ID in 1:4 is required.');
end
names = {'l6','l6','inhibition','inhibition'};
signs = [-1 1 -1 1];
stem = sprintf('%s_sign_%+d_result.mat',names{taskId},signs(taskId));
stem = strrep(stem,'+','p');
stem = strrep(stem,'-','m');
setupLoaded = load(setupFile,'setup');
branchLoaded = load(fullfile(coarseRoot,stem),'result');
result = l6ns_global_bifurcation_refine(setupLoaded.setup, ...
    branchLoaded.result,outputDir); %#ok<NASGU>
