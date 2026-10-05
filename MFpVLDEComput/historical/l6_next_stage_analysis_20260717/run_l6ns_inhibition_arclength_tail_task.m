% Refine the inhibition continuation tail after rejecting a branch-switch jump.

setupFile = getenv('L6NS_GLOBAL_BIFURCATION_SETUP');
sourceRoot = getenv('L6NS_GLOBAL_BIFURCATION_ARCLENGTH_ROOT');
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_FINAL_BRANCH_ROOT');
if isempty(setupFile) || isempty(sourceRoot) || isempty(outputDir)
    error('Setup, source arclength, and final branch roots are required.');
end
taskId = str2double(getenv('L6NS_BRANCH_TASK_ID'));
if ~ismember(taskId,[3 4])
    error('The inhibition tail task ID must be 3 or 4.');
end
signs = [-1 1];
branchSign = signs(taskId-2);
stem = sprintf('inhibition_sign_%+d_result.mat',branchSign);
stem = strrep(stem,'+','p');
stem = strrep(stem,'-','m');
setupLoaded = load(setupFile,'setup');
branchLoaded = load(fullfile(sourceRoot,stem),'result');
source = branchLoaded.result;

distance = source.Table.branchDistance;
jump = find(distance(2:end)<0.2*distance(1:end-1) & ...
    distance(1:end-1)>1,1,'first');
if isempty(jump)
    error('No persistent-branch jump was detected in task %d.',taskId);
end
source.Table = source.Table(1:jump,:);
source.States = source.States(1:jump);
options = struct('MinimumStep',5e-4,'MaximumStep',0.12, ...
    'MaximumPoints',120,'ReturnDistanceRatio',0.002);
result = l6ns_global_bifurcation_arclength(setupLoaded.setup,source, ...
    outputDir,options); %#ok<NASGU>
