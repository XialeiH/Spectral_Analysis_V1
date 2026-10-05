function run_l6_ei_grid_task()
% Evaluate one independent (beta6,betaEI) moved-fixed-point condition.

localTaskId=str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isfinite(localTaskId)
    taskOffset=str2double(getenv('L6EI_TASK_OFFSET'));
    if ~isfinite(taskOffset); taskOffset=0; end
    taskListIndex=localTaskId+taskOffset;
    taskListFile=getenv('L6EI_TASK_LIST_FILE');
    if ~isempty(taskListFile)
        taskIds=readmatrix(taskListFile,'FileType','text');
        if taskListIndex<1 || taskListIndex>numel(taskIds)
            error('L6EI:TaskListIndex', ...
                'Missing-task list index must be in 1:%d.',numel(taskIds));
        end
        taskId=taskIds(taskListIndex);
    else
        taskId=taskListIndex;
    end
else
    taskId=str2double(getenv('L6EI_TASK_ID'));
end
gridStep=local_env_number('L6EI_GRID_STEP',0.005);
beta6Maximum=local_env_number('L6EI_BETA6_MAXIMUM',0.35);
betaEIMaximum=local_env_number('L6EI_BETAEI_MAXIMUM',0.2);
beta6Count=round(beta6Maximum/gridStep)+1;
betaEICount=round(betaEIMaximum/gridStep)+1;
totalCount=beta6Count*betaEICount;
if ~isfinite(taskId) || taskId<1 || taskId>totalCount || ...
        taskId~=round(taskId)
    error('L6EI:TaskId','Task id must be an integer in 1:%d.',totalCount);
end

setupFile=getenv('L6EI_SETUP_FILE');
outputRoot=getenv('L6EI_OUTPUT_ROOT');
if ~isfile(setupFile) || isempty(outputRoot)
    error('L6EI:Environment','Setup file and output root are required.');
end
loaded=load(setupFile,'setup');
context=loaded.setup.Context;
baseline=context.FixedPoint(:);
if abs(context.CWeight-0.3077)>1e-12
    error('L6EI:HCWeight', ...
        'Context C weight does not match canonical HC_norm_diff.');
end

beta6Index=floor((taskId-1)/betaEICount)+1;
betaEIIndex=mod(taskId-1,betaEICount)+1;
beta6=(beta6Index-1)*gridStep;
betaEI=(betaEIIndex-1)*gridStep;
gain6=1+beta6;
gainEI=1+betaEI;

phi=@(state)l6ns_phi_l6_ei(state,gain6,gainEI,context);
fixed=real_tuning_fixed_point(phi,baseline,context.RelaxationP);
finiteEndpoint=all(isfinite(fixed.State)) && isfinite(fixed.Residual);
accepted=finiteEndpoint && fixed.Residual<=1e-5 && ...
    max(abs(fixed.State))<fixed.RateLimitHz;
n=numel(baseline)/3;
hcNorm=HC_norm_diff( ...
    fixed.State(1:n),fixed.State(n+(1:n)),fixed.State(2*n+(1:n)), ...
    baseline(1:n),baseline(n+(1:n)),baseline(2*n+(1:n)), ...
    0.3077,0.8,0.2);
minimumRate=min(fixed.State);
maximumRate=max(fixed.State);

pointsRoot=fullfile(outputRoot,'points');
if ~exist(pointsRoot,'dir'); mkdir(pointsRoot); end
outputFile=fullfile(pointsRoot,sprintf('point_%05d.tsv',taskId));
temporaryFile=[outputFile '.tmp'];
fileId=fopen(temporaryFile,'w');
if fileId<0; error('L6EI:Output','Cannot open %s.',temporaryFile); end
cleanup=onCleanup(@()fclose(fileId));
fprintf(fileId,[ ...
    '%d\t%d\t%d\t%.15g\t%.15g\t%.15g\t%.15g\t%d\t%d\t' ...
    '%.15g\t%d\t%.15g\t%.15g\t%.15g\n'], ...
    taskId,beta6Index,betaEIIndex,beta6,betaEI,gain6,gainEI, ...
    fixed.Converged,accepted,fixed.Residual,fixed.Iterations,hcNorm, ...
    minimumRate,maximumRate);
clear cleanup
movefile(temporaryFile,outputFile,'f');
fprintf(['L6-EI grid %d/%d: beta6=%.3f betaEI=%.3f strict=%d ' ...
    'accepted=%d HC=%.9g residual=%.3g iterations=%d.\n'], ...
    taskId,totalCount,beta6,betaEI,fixed.Converged,accepted,hcNorm, ...
    fixed.Residual,fixed.Iterations);
end

function value=local_env_number(name,defaultValue)
value=str2double(getenv(name));
if ~isfinite(value); value=defaultValue; end
end
