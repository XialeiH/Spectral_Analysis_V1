function run_pair_grid_ode_fallback_task()
% Re-evaluate one nonconverged signed L6-excitation point with ode45.

fallbackIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(fallbackIndex)
    fallbackIndex = str2double(getenv('PAIR_GRID_ODE_INDEX'));
end
taskListFile = getenv('PAIR_GRID_ODE_TASK_LIST');
taskIds = readmatrix(taskListFile,'FileType','text');
taskIds = taskIds(isfinite(taskIds));
if ~isfinite(fallbackIndex) || fallbackIndex<1 || ...
        fallbackIndex>numel(taskIds) || fallbackIndex~=round(fallbackIndex)
    error('PairGridODE:Index','Fallback index must be in 1:%d.',numel(taskIds));
end
taskId = taskIds(fallbackIndex);

setupFile = getenv('PAIR_GRID_SETUP_FILE');
outputRoot = getenv('PAIR_GRID_ODE_OUTPUT_ROOT');
loaded = load(setupFile,'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
gridStep = 0.005;
gridMaximum = 0.4;
gridCount = round(gridMaximum/gridStep)+1;
betaXIndex = floor((taskId-1)/gridCount)+1;
betaYIndex = mod(taskId-1,gridCount)+1;
beta6 = (betaXIndex-1)*gridStep;
betaE = -gridMaximum+(betaYIndex-1)*gridStep;

phi = @(state)l6ns_phi_three_pathway( ...
    state,1+beta6,1+betaE,1,context);
fixed = ode_tuning_fixed_point(phi,baseline);
n = numel(baseline)/3;
hcNorm = HC_norm_diff( ...
    fixed.State(1:n),fixed.State(n+(1:n)),fixed.State(2*n+(1:n)), ...
    baseline(1:n),baseline(n+(1:n)),baseline(2*n+(1:n)), ...
    0.3077,0.8,0.2);

pointsRoot = fullfile(outputRoot,'points');
if ~exist(pointsRoot,'dir'); mkdir(pointsRoot); end
outputFile = fullfile(pointsRoot,sprintf('ode_point_%05d.tsv',taskId));
temporaryFile = [outputFile '.tmp'];
fileId = fopen(temporaryFile,'w');
if fileId<0
    error('PairGridODE:Output','Cannot open %s.',temporaryFile);
end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,['%d\t%d\t%d\t%.15g\t%.15g\t%d\t%.15g\t%.15g\t' ...
    '%.15g\t%.15g\t%.15g\t%.15g\t%s\n'],taskId,betaXIndex,betaYIndex,beta6,betaE, ...
    fixed.Converged,fixed.Residual,fixed.IntrinsicTime,fixed.AcceptedSteps, ...
    hcNorm,min(fixed.State),max(fixed.State),char(fixed.Termination));
clear cleanup
movefile(temporaryFile,outputFile,'f');
fprintf(['ODE fallback %d/%d, task %d: beta6=%.3f betaE=%.3f, ' ...
    'converged=%d residual=%.3g intrinsicTime=%.3g steps=%d HC=%.9g.\n'], ...
    fallbackIndex,numel(taskIds),taskId,beta6,betaE,fixed.Converged, ...
    fixed.Residual,fixed.IntrinsicTime,fixed.AcceptedSteps,hcNorm);
end
