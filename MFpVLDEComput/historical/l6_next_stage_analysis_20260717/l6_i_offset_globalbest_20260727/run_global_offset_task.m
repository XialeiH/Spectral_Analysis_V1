function run_global_offset_task()
% Find the global best scalar betaI for one fixed beta6.

taskId=str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId); taskId=str2double(getenv('GLOBAL_OFFSET_TASK_ID')); end
smoke=strcmp(getenv('GLOBAL_OFFSET_SMOKE'),'1');
cases=global_offset_case_table(smoke);
if ~isfinite(taskId) || taskId<1 || taskId>height(cases)
    error('GlobalOffset:TaskId','Task id must be in 1:%d.',height(cases));
end
spec=cases(round(taskId),:);
setupFile=getenv('GLOBAL_OFFSET_SETUP_FILE');
outputRoot=getenv('GLOBAL_OFFSET_OUTPUT_ROOT');
if ~isfile(setupFile) || isempty(outputRoot)
    error('GlobalOffset:Environment','Setup file and output root are required.');
end
loaded=load(setupFile,'setup');
context=loaded.setup.Context;
baseline=context.FixedPoint(:);
beta6=spec.beta6;
hcWC=0.3077;
hcAE=32/(32+8);
hcAI=8/(32+8);
if abs(context.CWeight-hcWC)>1e-12
    error('GlobalOffset:HCWeight', ...
        'Context C weight %.17g differs from canonical HC weight %.17g.', ...
        context.CWeight,hcWC);
end

cacheBetaI=zeros(0,1);
cacheHC=zeros(0,1);
cacheConverged=false(0,1);
cacheResidual=zeros(0,1);
cacheIterations=zeros(0,1);
cacheMinimumRate=zeros(0,1);
cacheMaximumRate=zeros(0,1);
cacheStage=strings(0,1);
activeStage="coarse";

% Search betaI directly on one absolute interval shared by every beta6.
% Ratios are computed only after all independently optimized points exist.
betaILower=-0.25;
betaIUpper=0.25;
coarseCount=25;
maximumExpansions=3;
expansionCount=0;
while true
    activeStage="coarse";
    coarseBetaI=linspace(betaILower,betaIUpper,coarseCount)';
    coarseHC=nan(size(coarseBetaI));
    for index=1:numel(coarseBetaI)
        coarseHC(index)=objective(coarseBetaI(index));
    end
    [~,coarseBestIndex]=min(coarseHC);
    if coarseBestIndex>1 && coarseBestIndex<numel(coarseBetaI)
        break
    end
    if expansionCount>=maximumExpansions
        error('GlobalOffset:Boundary', ...
            'Best coarse point remains on a boundary after %d expansions.', ...
            maximumExpansions);
    end
    width=betaIUpper-betaILower;
    if coarseBestIndex==1
        betaILower=max(-0.99,betaILower-width);
    else
        betaIUpper=betaIUpper+width;
    end
    expansionCount=expansionCount+1;
end

localIndices=find(coarseHC(2:end-1)<=coarseHC(1:end-2) & ...
    coarseHC(2:end-1)<=coarseHC(3:end) & coarseHC(2:end-1)<1e5)+1;
localIndices=unique([localIndices;coarseBestIndex]);
localMinimumCount=numel(localIndices);
activeStage="refine";
refinedBetaI=zeros(localMinimumCount,1);
for index=1:localMinimumCount
    location=localIndices(index);
    lower=coarseBetaI(location-1);
    upper=coarseBetaI(location+1);
    options=optimset('Display','off','TolX',1e-10, ...
        'MaxIter',100,'MaxFunEvals',160);
    refinedBetaI(index)=fminbnd(@objective,lower,upper,options);
    objective(refinedBetaI(index));
end
candidates=unique([coarseBetaI;refinedBetaI;0;0.515*beta6]);
candidateHC=nan(size(candidates));
activeStage="candidate";
for index=1:numel(candidates)
    candidateHC(index)=objective(candidates(index));
end
[~,bestIndex]=min(candidateHC);
betaIBest=candidates(bestIndex);
boundaryGap=min(coarseHC([1 end]))-candidateHC(bestIndex);
boundaryMinimum=abs(betaIBest-coarseBetaI(1))<1e-10 || ...
    abs(betaIBest-coarseBetaI(end))<1e-10;

verified=evaluate_pair(beta6,betaIBest,baseline,context);
if ~verified.Converged
    error('GlobalOffset:Verification','Independent fixed-point verification failed.');
end
caseTag=sprintf('%03d_beta6_%s',taskId,number_tag(beta6));
caseDir=fullfile(outputRoot,caseTag);
if ~exist(caseDir,'dir'); mkdir(caseDir); end

[sortedBetaI,order]=sort(cacheBetaI);
profile=table(sortedBetaI,cacheHC(order),cacheConverged(order), ...
    cacheResidual(order),cacheIterations(order),cacheMinimumRate(order), ...
    cacheMaximumRate(order),cacheStage(order), ...
    'VariableNames',{'betaI','HCnorm','converged','fixedPointResidual', ...
    'iterations','minimumRateHz','maximumRateHz','searchStage'});
writetable(profile,fullfile(caseDir,'global_search_profile.tsv'), ...
    'FileType','text','Delimiter','\t');

summary=table(taskId,beta6,betaIBest,1+beta6,1+betaIBest, ...
    verified.Converged,string(verified.Termination),verified.Iterations, ...
    verified.Residual,verified.HCnorm,verified.RelativeShift, ...
    verified.MinimumRateHz,verified.MaximumRateHz,height(profile), ...
    betaILower,betaIUpper,expansionCount,localMinimumCount,boundaryGap, ...
    boundaryMinimum,hcWC,hcAE,hcAI, ...
    'VariableNames',{'taskId','beta6','betaI','gain6','gainI','converged', ...
    'termination','iterations','fixedPointResidual','HCnorm', ...
    'relativeFixedPointShift','minimumRateHz','maximumRateHz', ...
    'searchEvaluations','betaILower','betaIUpper','expansionCount', ...
    'localMinimumCount','boundaryGap','boundaryMinimum', ...
    'HCwC','HCaE','HCaI'});
writetable(summary,fullfile(caseDir,'case_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
fixedPoint=verified.State;
save(fullfile(caseDir,'case_result.mat'),'summary','profile','fixedPoint','-v7.3');
fprintf(['Global offset %d: beta6=%+.5f betaI=%+.10f HC=%.9g ' ...
    'eval=%d betaIRange=[%.3g %.3g] boundaryGap=%.6g.\n'], ...
    taskId,beta6,betaIBest,verified.HCnorm,height(profile), ...
    betaILower,betaIUpper,boundaryGap);

    function value=objective(betaI)
        match=find(abs(cacheBetaI-betaI)<1e-12,1);
        if ~isempty(match)
            value=cacheHC(match);
            return
        end
        one=evaluate_pair(beta6,betaI,baseline,context);
        value=one.HCnorm;
        if ~one.Converged || ~isfinite(value); value=1e6; end
        cacheBetaI(end+1,1)=betaI;
        cacheHC(end+1,1)=value;
        cacheConverged(end+1,1)=one.Converged;
        cacheResidual(end+1,1)=one.Residual;
        cacheIterations(end+1,1)=one.Iterations;
        cacheMinimumRate(end+1,1)=one.MinimumRateHz;
        cacheMaximumRate(end+1,1)=one.MaximumRateHz;
        cacheStage(end+1,1)=activeStage;
    end
end

function result=evaluate_pair(beta6,betaI,baseline,context)
gain6=1+beta6;
gainI=1+betaI;
if gainI<=0
    result=struct('Converged',false,'Termination',"nonpositive_gain", ...
        'Iterations',0,'Residual',Inf,'State',baseline,'MinimumRateHz',NaN, ...
        'MaximumRateHz',NaN,'HCnorm',Inf,'RelativeShift',Inf);
    return
end
phi=@(state)l6ns_phi(state,1-gain6,context,[1 1],gainI,'true');
fixed=real_tuning_fixed_point(phi,baseline,context.RelaxationP);
n=numel(baseline)/3;
result=fixed;
result.HCnorm=HC_norm_diff( ...
    fixed.State(1:n),fixed.State(n+(1:n)),fixed.State(2*n+(1:n)), ...
    baseline(1:n),baseline(n+(1:n)),baseline(2*n+(1:n)), ...
    0.3077,32/(32+8),8/(32+8));
result.RelativeShift=norm(fixed.State-baseline)/max(norm(baseline),eps);
end

function tag=number_tag(value)
tag=sprintf('%+.5f',value);
tag=strrep(tag,'+','p'); tag=strrep(tag,'-','m'); tag=strrep(tag,'.','p');
end
