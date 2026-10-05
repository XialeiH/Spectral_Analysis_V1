function run_pair_grid_task()
% Evaluate one independent moved-fixed-point pathway-pair grid condition.

localTaskId=str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isfinite(localTaskId)
    taskOffset=str2double(getenv('PAIR_GRID_TASK_OFFSET'));
    if ~isfinite(taskOffset); taskOffset=0; end
    taskId=localTaskId+taskOffset;
else
    taskId=str2double(getenv('PAIR_GRID_TASK_ID'));
end
gridStep=0.005;
gridMaximum=0.4;
gridCount=round(gridMaximum/gridStep)+1;
totalCount=gridCount^2;
if ~isfinite(taskId) || taskId<1 || taskId>totalCount || ...
        taskId~=round(taskId)
    error('PairGrid:TaskId','Task id must be an integer in 1:%d.',totalCount);
end

setupFile=getenv('PAIR_GRID_SETUP_FILE');
outputRoot=getenv('PAIR_GRID_OUTPUT_ROOT');
pairType=string(getenv('PAIR_GRID_PAIR'));
if ~isfile(setupFile) || isempty(outputRoot)
    error('PairGrid:Environment','Setup file and output root are required.');
end
if ~any(pairType==["l6_excitation","excitation_inhibition"])
    error('PairGrid:PairType','Unsupported pair type: %s.',pairType);
end
loaded=load(setupFile,'setup');
context=loaded.setup.Context;
baseline=context.FixedPoint(:);
if abs(context.CWeight-0.3077)>1e-12
    error('PairGrid:HCWeight', ...
        'Context C weight does not match canonical HC_norm_diff.');
end

betaXIndex=floor((taskId-1)/gridCount)+1;
betaYIndex=mod(taskId-1,gridCount)+1;
betaX=(betaXIndex-1)*gridStep;
beta6=0;
betaE=0;
betaI=0;
switch pairType
    case "l6_excitation"
        betaY=-gridMaximum+(betaYIndex-1)*gridStep;
        beta6=betaX;
        betaE=betaY;
    case "excitation_inhibition"
        betaY=(betaYIndex-1)*gridStep;
        betaE=betaX;
        betaI=betaY;
end
gain6=1+beta6;
gainE=1+betaE;
gainI=1+betaI;

phi=@(state)l6ns_phi_three_pathway(state,gain6,gainE,gainI,context);
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
if fileId<0
    error('PairGrid:Output','Cannot open %s.',temporaryFile);
end
cleanup=onCleanup(@()fclose(fileId));
fprintf(fileId,['%d\t%d\t%d\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t' ...
    '%d\t%.15g\t%d\t%.15g\t%.15g\t%.15g\n'], ...
    taskId,betaXIndex,betaYIndex,betaX,betaY,beta6,betaE,betaI, ...
    fixed.Converged,fixed.Residual,fixed.Iterations,hcNorm, ...
    minimumRate,maximumRate);
clear cleanup
movefile(temporaryFile,outputFile,'f');
fprintf(['Pair grid %s %d/%d: beta6=%.3f betaE=%.3f betaI=%.3f ' ...
    'converged=%d HC=%.9g residual=%.3g iterations=%d.\n'], ...
    pairType,taskId,totalCount,beta6,betaE,betaI,fixed.Converged, ...
    hcNorm,fixed.Residual,fixed.Iterations);
end
