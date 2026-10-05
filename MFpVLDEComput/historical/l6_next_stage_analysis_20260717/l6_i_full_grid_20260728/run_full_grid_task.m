function run_full_grid_task()
% Evaluate one independent (beta6,betaI) moved-fixed-point grid condition.

localTaskId=str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isfinite(localTaskId)
    taskOffset=str2double(getenv('FULL_GRID_TASK_OFFSET'));
    if ~isfinite(taskOffset); taskOffset=0; end
    taskId=localTaskId+taskOffset;
else
    taskId=str2double(getenv('FULL_GRID_TASK_ID'));
end
gridStep=0.002;
gridMaximum=0.4;
gridCount=round(gridMaximum/gridStep)+1;
totalCount=gridCount^2;
if ~isfinite(taskId) || taskId<1 || taskId>totalCount || taskId~=round(taskId)
    error('FullGrid:TaskId','Task id must be an integer in 1:%d.',totalCount);
end

setupFile=getenv('FULL_GRID_SETUP_FILE');
outputRoot=getenv('FULL_GRID_OUTPUT_ROOT');
if ~isfile(setupFile) || isempty(outputRoot)
    error('FullGrid:Environment','Setup file and output root are required.');
end
loaded=load(setupFile,'setup');
context=loaded.setup.Context;
baseline=context.FixedPoint(:);
if abs(context.CWeight-0.3077)>1e-12
    error('FullGrid:HCWeight','Context C weight does not match canonical HC_norm_diff.');
end

beta6Index=floor((taskId-1)/gridCount)+1;
betaIIndex=mod(taskId-1,gridCount)+1;
beta6=(beta6Index-1)*gridStep;
betaI=(betaIIndex-1)*gridStep;
gain6=1+beta6;
gainI=1+betaI;

phi=@(state)l6ns_phi(state,1-gain6,context,[1 1],gainI,'true');
fixed=real_tuning_fixed_point(phi,baseline,context.RelaxationP);
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
if fileId<0; error('FullGrid:Output','Cannot open %s.',temporaryFile); end
cleanup=onCleanup(@()fclose(fileId));
fprintf(fileId,['%d\t%d\t%d\t%.15g\t%.15g\t%d\t%.15g\t%d\t' ...
    '%.15g\t%.15g\t%.15g\n'], ...
    taskId,beta6Index,betaIIndex,beta6,betaI,fixed.Converged, ...
    fixed.Residual,fixed.Iterations,hcNorm,minimumRate,maximumRate);
clear cleanup
movefile(temporaryFile,outputFile,'f');
fprintf(['Full grid %d/%d: beta6=%.3f betaI=%.3f converged=%d ' ...
    'HC=%.9g residual=%.3g iterations=%d.\n'], ...
    taskId,totalCount,beta6,betaI,fixed.Converged,hcNorm, ...
    fixed.Residual,fixed.Iterations);
end
