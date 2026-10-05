function run_l6_i_offset_task()
% Optimize beta_I for one beta_6 using real tuning and regular iteration.

taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId); taskId=str2double(getenv('OFFSET_TASK_ID')); end
smoke = strcmp(getenv('OFFSET_SMOKE'),'1');
cases = offset_case_table(smoke);
if ~isfinite(taskId) || taskId<1 || taskId>height(cases)
    error('Offset:TaskId','Task id must be in 1:%d.',height(cases));
end
spec = cases(round(taskId),:);
setupFile = getenv('OFFSET_SETUP_FILE');
outputRoot = getenv('OFFSET_OUTPUT_ROOT');
if ~isfile(setupFile) || isempty(outputRoot)
    error('Offset:Environment','OFFSET_SETUP_FILE and OFFSET_OUTPUT_ROOT are required.');
end
loaded = load(setupFile,'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
beta6 = spec.beta6;

cacheBetaI = zeros(0,1);
cacheHC = zeros(0,1);
cacheConverged = false(0,1);
cacheResidual = zeros(0,1);
cacheIterations = zeros(0,1);
cacheMinimumRate = zeros(0,1);
cacheMaximumRate = zeros(0,1);

if beta6==0
    betaICandidate = 0;
else
    lowerBound = 0;
    upperBound = max(0.005,1.25*beta6);
    options = optimset('Display','off','TolX',1e-7, ...
        'MaxIter',40,'MaxFunEvals',50);
    betaIFmin = fminbnd(@objective,lowerBound,upperBound,options);
    candidates = unique([lowerBound; betaIFmin; upperBound; 0.515*beta6]);
    candidateHC = nan(size(candidates));
    for index=1:numel(candidates)
        candidateHC(index)=objective(candidates(index));
    end
    [~,bestIndex] = min(candidateHC);
    betaICandidate = candidates(bestIndex);
end

% Independent verification: recompute from the baseline without using cache.
verified = evaluate_pair(beta6,betaICandidate,baseline,context);
caseTag = sprintf('%03d_beta6_%s',taskId,number_tag(beta6));
caseDir = fullfile(outputRoot,caseTag);
if ~exist(caseDir,'dir'); mkdir(caseDir); end

trace = table(cacheBetaI,cacheHC,cacheConverged,cacheResidual, ...
    cacheIterations,cacheMinimumRate,cacheMaximumRate, ...
    'VariableNames',{'betaI','HCnorm','converged','fixedPointResidual', ...
    'iterations','minimumRateHz','maximumRateHz'});
writetable(trace,fullfile(caseDir,'optimization_trace.tsv'), ...
    'FileType','text','Delimiter','\t');

summary = table(taskId,beta6,betaICandidate,1+beta6,1+betaICandidate, ...
    verified.Converged,string(verified.Termination),verified.Iterations, ...
    verified.Residual,verified.HCnorm,verified.RelativeShift, ...
    verified.MinimumRateHz,verified.MaximumRateHz,height(trace), ...
    verified.HCnorm<0.3, ...
    'VariableNames',{'taskId','beta6','betaI','gain6','gainI','converged', ...
    'termination','iterations','fixedPointResidual','HCnorm', ...
    'relativeFixedPointShift','minimumRateHz','maximumRateHz', ...
    'optimizationEvaluations','isUnmoved'});
writetable(summary,fullfile(caseDir,'case_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
fixedPoint = verified.State;
save(fullfile(caseDir,'case_result.mat'),'summary','trace','fixedPoint','-v7.3');
fprintf(['Offset task %d: beta6=%+.5f betaI=%+.8f HC=%.6g ' ...
    'converged=%d residual=%.3e evaluations=%d.\n'],taskId,beta6, ...
    betaICandidate,verified.HCnorm,verified.Converged,verified.Residual,height(trace));

    function value = objective(betaI)
        match = find(abs(cacheBetaI-betaI)<1e-12,1);
        if ~isempty(match)
            value = cacheHC(match);
            return
        end
        one = evaluate_pair(beta6,betaI,baseline,context);
        value = one.HCnorm;
        if ~one.Converged || ~isfinite(value); value=1e6; end
        cacheBetaI(end+1,1)=betaI;
        cacheHC(end+1,1)=value;
        cacheConverged(end+1,1)=one.Converged;
        cacheResidual(end+1,1)=one.Residual;
        cacheIterations(end+1,1)=one.Iterations;
        cacheMinimumRate(end+1,1)=one.MinimumRateHz;
        cacheMaximumRate(end+1,1)=one.MaximumRateHz;
    end
end

function result = evaluate_pair(beta6,betaI,baseline,context)
gain6 = 1+beta6;
gainI = 1+betaI;
phi = @(state)l6ns_phi(state,1-gain6,context,[1 1],gainI,'true');
fixed = real_tuning_fixed_point(phi,baseline,context.RelaxationP);
n = numel(baseline)/3;
dS = fixed.State(1:n)-baseline(1:n);
dC = fixed.State(n+(1:n))-baseline(n+(1:n));
dI = fixed.State(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-context.CWeight)*dS+context.CWeight*dC;
hc = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
result = fixed;
result.HCnorm = hc;
result.RelativeShift = norm(fixed.State-baseline)/max(norm(baseline),eps);
end

function tag = number_tag(value)
tag=sprintf('%+.5f',value);
tag=strrep(tag,'+','p'); tag=strrep(tag,'-','m'); tag=strrep(tag,'.','p');
end
