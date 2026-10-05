function run_fpp_offset_task()
% Match the baseline leading ODE rate using simultaneous FPP L6 and I gains.

taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId); taskId=str2double(getenv('FPP_OFFSET_TASK_ID')); end
smoke = strcmp(getenv('FPP_OFFSET_SMOKE'),'1');
cases = fpp_offset_case_table(smoke);
if ~isfinite(taskId) || taskId<1 || taskId>height(cases)
    error('FppOffset:TaskId','Task id must be in 1:%d.',height(cases));
end
spec = cases(round(taskId),:);
setupFile = getenv('FPP_OFFSET_SETUP_FILE');
theoryFile = getenv('FPP_OFFSET_THEORY_FILE');
outputRoot = getenv('FPP_OFFSET_OUTPUT_ROOT');
if ~isfile(setupFile) || ~isfile(theoryFile) || isempty(outputRoot)
    error('FppOffset:Environment', ...
        'Setup, theory, and output environment variables are required.');
end

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
fixed = context.FixedPoint(:);
pathway = setup.Pathway;
cfg = setup.Config;
if pathway.JINonISourceFraction>1e-12 || pathway.J6ISourceFraction>1e-12
    error('FppOffset:PathwayDefinition', ...
        ['Expected J_I to contain only inhibitory-source columns and J_6 ' ...
        'to contain no inhibitory-source columns; measured %.3e and %.3e.'], ...
        pathway.JINonISourceFraction,pathway.J6ISourceFraction);
end
theory = readtable(theoryFile,'FileType','text','Delimiter','\t');
theorySlope = theory.constantStabilitySlope(1);
beta6 = spec.beta6;
gamma6 = 1+beta6;

[baselineAlpha,baselineLambda,baselineRight,~,baselineS6,baselineSI] = ...
    local_leading_mode(pathway.JBaseline,pathway,cfg);
if abs(baselineAlpha-theory.baselineLambdaReal(1))>2e-7
    error('FppOffset:BaselineAlpha', ...
        'Baseline alpha %.12g differs from theory %.12g.', ...
        baselineAlpha,theory.baselineLambdaReal(1));
end

betaI = theorySlope*beta6;
traceIteration = zeros(0,1);
traceBetaI = zeros(0,1);
traceAlpha = zeros(0,1);
traceError = zeros(0,1);
traceS6 = zeros(0,1);
traceSI = zeros(0,1);
maximumNewtonIterations = 8;
alphaTolerance = 2e-9;
for iteration = 1:maximumNewtonIterations
    gammaI = 1+betaI;
    jacobian = pathway.JRest+gamma6*pathway.J6+gammaI*pathway.JI;
    [alpha,lambda,~,~,s6,sI] = local_leading_mode(jacobian,pathway,cfg);
    alphaError = alpha-baselineAlpha;
    traceIteration(end+1,1)=iteration;
    traceBetaI(end+1,1)=betaI;
    traceAlpha(end+1,1)=alpha;
    traceError(end+1,1)=alphaError;
    traceS6(end+1,1)=s6;
    traceSI(end+1,1)=sI;
    if abs(alphaError)<=alphaTolerance; break; end
    if ~isfinite(sI) || abs(sI)<1e-8
        error('FppOffset:Derivative','Inhibitory eigenvalue sensitivity is singular.');
    end
    candidate = betaI-alphaError/sI;
    if ~isfinite(candidate) || candidate<0
        candidate = max(0,0.5*betaI);
    end
    betaI = candidate;
end
gammaI = 1+betaI;
jacobian = pathway.JRest+gamma6*pathway.J6+gammaI*pathway.JI;
[alpha,lambda,~,~,s6,sI] = local_leading_mode(jacobian,pathway,cfg);
alphaError = alpha-baselineAlpha;
converged = abs(alphaError)<=1e-7;
if ~converged
    error('FppOffset:NoConvergence', ...
        'Constant-stability solve ended with alpha error %.3e.',alphaError);
end

phi = @(state)l6ns_phi(state,1-gamma6,context,[1 1],gammaI,'fpp');
fixedPointResidual = norm(phi(fixed)-fixed)/max(norm(fixed),eps);
fixedPointHCnorm = local_hcnorm(fixed,fixed,context.CWeight);

rng(7000+taskId);
direction = randn(size(fixed)); direction=direction/norm(direction);
step = 2e-6*max(1,norm(fixed)/sqrt(numel(fixed)));
numeric = (phi(fixed+step*direction)-phi(fixed-step*direction))/(2*step);
analyticJvpError = norm(numeric-jacobian*direction)/max(norm(jacobian*direction),eps);

[odeTrace,odeSummary] = local_run_ode(phi,fixed,baselineRight,context, ...
    cfg.TauMs,beta6,betaI,"matched");
odeSummary.PredictedDecayPerMs = (alpha-1)/cfg.TauMs;
caseTag = sprintf('%03d_beta6_%s',taskId,local_number_tag(beta6));
caseDir = fullfile(outputRoot,caseTag);
if ~exist(caseDir,'dir'); mkdir(caseDir); end
writetable(odeTrace,fullfile(caseDir,'ode_trace.tsv'), ...
    'FileType','text','Delimiter','\t');

if abs(beta6-max(cases.beta6))<1e-12 && beta6>0
    phiL6Only = @(state)l6ns_phi(state,1-gamma6,context,[1 1],1,'fpp');
    [l6OnlyTrace,l6OnlySummary] = local_run_ode(phiL6Only,fixed, ...
        baselineRight,context,cfg.TauMs,beta6,0,"l6_only");
    writetable(l6OnlyTrace,fullfile(caseDir,'ode_trace_l6_only.tsv'), ...
        'FileType','text','Delimiter','\t');
else
    l6OnlySummary = struct('FittedDecayPerMs',NaN,'PredictedDecayPerMs',NaN, ...
        'FinalHCnorm',NaN);
end

trace = table(traceIteration,traceBetaI,traceAlpha,traceError,traceS6,traceSI, ...
    'VariableNames',{'iteration','betaI','leadingAlpha','alphaError','s6','sI'});
writetable(trace,fullfile(caseDir,'newton_trace.tsv'), ...
    'FileType','text','Delimiter','\t');

summary = table(taskId,beta6,betaI,gamma6,gammaI,theorySlope, ...
    betaI/max(beta6,eps),baselineAlpha,real(baselineLambda), ...
    imag(baselineLambda),alpha,real(lambda),imag(lambda),alphaError, ...
    converged,height(trace),s6,sI,-s6/sI,baselineS6,baselineSI, ...
    fixedPointResidual,fixedPointHCnorm,analyticJvpError,cfg.TauMs, ...
    odeSummary.PredictedDecayPerMs,odeSummary.FittedDecayPerMs, ...
    odeSummary.FittedDecayPerMs-odeSummary.PredictedDecayPerMs, ...
    odeSummary.FinalHCnorm,l6OnlySummary.FittedDecayPerMs, ...
    l6OnlySummary.FinalHCnorm,pathway.JINonISourceFraction, ...
    pathway.J6ISourceFraction, ...
    'VariableNames',{'taskId','beta6','betaI','gamma6','gammaI', ...
    'theoreticalSlope','betaIRatio','baselineMaxRealLambda', ...
    'baselineLeadingReal','baselineLeadingImag','matchedMaxRealLambda', ...
    'matchedLeadingReal','matchedLeadingImag','matchingError','converged', ...
    'newtonIterations','matchedS6','matchedSI','matchedLocalSlope', ...
    'baselineS6','baselineSI','fixedPointResidual','fixedPointHCnorm', ...
    'analyticJvpRelativeError','tauMs','predictedDecayPerMs', ...
    'fittedDecayPerMs','fittedMinusPredictedDecayPerMs','finalHCnorm', ...
    'l6OnlyFittedDecayPerMs','l6OnlyFinalHCnorm', ...
    'JINonISourceFraction','J6ISourceFraction'});
writetable(summary,fullfile(caseDir,'case_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(caseDir,'case_result.mat'),'summary','trace','odeTrace','-v7.3');
fprintf(['FPP offset %d: beta6=%+.5f betaI=%+.8f ratio=%.6f ' ...
    'alphaError=%+.3e fixedResidual=%.3e ODErateError=%+.3e.\n'], ...
    taskId,beta6,betaI,betaI/max(beta6,eps),alphaError, ...
    fixedPointResidual,summary.fittedMinusPredictedDecayPerMs);
end

function [alpha,lambda,r,l,s6,sI] = local_leading_mode(jacobian,pathway,cfg)
count = 6;
options = struct('tol',cfg.EigsTolerance,'maxit',cfg.EigsMaxIterations, ...
    'p',min(max(cfg.EigsSubspaceDimension,2*count+8),size(jacobian,1)), ...
    'disp',0,'isreal',true);
try
    [rightPool,rightValues] = eigs(jacobian,count,'largestreal',options);
catch
    [rightPool,rightValues] = eigs(jacobian,count,'lr',options);
end
rightLambda = diag(rightValues);
[alpha,index] = max(real(rightLambda));
lambda = rightLambda(index);
r = rightPool(:,index); r=r/norm(r);
try
    [leftPool,leftValues] = eigs(jacobian',count,'largestreal',options);
catch
    [leftPool,leftValues] = eigs(jacobian',count,'lr',options);
end
leftLambda = diag(leftValues);
[~,leftIndex] = min(abs(leftLambda-conj(lambda)));
l = leftPool(:,leftIndex);
overlap = l'*r;
if abs(overlap)<1e-12
    error('FppOffset:Eigenvectors','Leading left/right overlap is singular.');
end
l = l/conj(overlap);
s6 = real(l'*(pathway.J6*r));
sI = real(l'*(pathway.JI*r));
end

function [trace,summary] = local_run_ode(phi,fixed,direction,context,tauMs, ...
        beta6,betaI,kind)
direction = real(direction(:)); direction=direction/norm(direction);
[~,index] = max(abs(direction));
if direction(index)<0; direction=-direction; end
unitHC = local_hcnorm(fixed+direction,fixed,context.CWeight);
amplitude = 0.01/max(unitHC,eps);
negative = direction<0;
if any(negative)
    amplitude=min(amplitude,0.5*min(fixed(negative)./(-direction(negative))));
end
initial = fixed+amplitude*direction;
timeMs = (0:5:500)';
rhs = @(~,state)(phi(state)-state)/tauMs;
options = odeset('RelTol',2e-6,'AbsTol',1e-8,'MaxStep',5);
solution = ode45(rhs,[timeMs(1) timeMs(end)],initial,options);
states = deval(solution,timeMs);
hc = nan(size(timeMs));
for index=1:numel(timeMs)
    hc(index)=local_hcnorm(states(:,index),fixed,context.CWeight);
end
fitMask = timeMs>=150 & timeMs<=400 & hc>1e-10;
fit = polyfit(timeMs(fitMask),log(hc(fitMask)),1);
fittedDecay = fit(1);
trace = table(repmat(string(kind),numel(timeMs),1), ...
    repmat(beta6,numel(timeMs),1),repmat(betaI,numel(timeMs),1), ...
    timeMs,hc,'VariableNames',{'kind','beta6','betaI','timeMs','HCnorm'});
summary = struct('FittedDecayPerMs',fittedDecay, ...
    'PredictedDecayPerMs',NaN,'FinalHCnorm',hc(end));
end

function value = local_hcnorm(state,baseline,wC)
n=numel(baseline)/3;
dS=state(1:n)-baseline(1:n);
dC=state(n+(1:n))-baseline(n+(1:n));
dI=state(2*n+(1:n))-baseline(2*n+(1:n));
dE=(1-wC)*dS+wC*dC;
value=sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end

function tag = local_number_tag(value)
tag=sprintf('%+.5f',value);
tag=strrep(tag,'+','p'); tag=strrep(tag,'-','m'); tag=strrep(tag,'.','p');
end
