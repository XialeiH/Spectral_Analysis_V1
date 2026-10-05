function run_real_tuning_task()
% Run one whole-pathway beta case with regular iteration and moved Jacobian.

taskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId)
    taskId = str2double(getenv('REAL_TUNING_TASK_ID'));
end
smoke = strcmp(getenv('REAL_TUNING_SMOKE'),'1');
cases = real_tuning_case_table(smoke);
if ~isfinite(taskId) || taskId<1 || taskId>height(cases)
    error('RealTuning:TaskId','Task id must be in 1:%d.',height(cases));
end
spec = cases(round(taskId),:);
setupFile = getenv('REAL_TUNING_SETUP_FILE');
outputRoot = getenv('REAL_TUNING_OUTPUT_ROOT');
if ~isfile(setupFile) || isempty(outputRoot)
    error('RealTuning:Environment','Setup file and output root are required.');
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
baseline = context.FixedPoint(:);
gain6 = 1;
gainI = 1;
if spec.pathway == "L6"
    gain6 = spec.pathwayGain;
else
    gainI = spec.pathwayGain;
end
phi = @(state)l6ns_phi(state,1-gain6,context,[1 1],gainI,'true');

fprintf('Real tuning task %d/%d: %s beta=%+.4f, gain=%.4f.\n', ...
    taskId,height(cases),spec.pathway,spec.beta,spec.pathwayGain);
fixedResult = real_tuning_fixed_point(phi,baseline,context.RelaxationP);

maxRealTrue = NaN;
leadingImagTrue = NaN;
maxRealFpp = NaN;
leadingImagFpp = NaN;
jvpRelativeError = NaN;
baselineJacobianRelativeError = NaN;
eigenvaluesTrue = [];
if fixedResult.Converged
    moved = fixedResult.State;
    J = real_tuning_true_jacobian(moved,context,gain6,gainI);
    [maxRealTrue,leadingImagTrue] = local_leading(J);
    jvpRelativeError = local_jvp_check(J,phi,moved,taskId);
    if spec.pathway == "L6"
        jFpp = setup.Pathway.JRest+gain6*setup.Pathway.J6+setup.Pathway.JI;
    else
        jFpp = setup.Pathway.JRest+setup.Pathway.J6+gainI*setup.Pathway.JI;
    end
    [maxRealFpp,leadingImagFpp] = local_leading(jFpp);
    if abs(spec.beta)<1e-12
        baselineJacobianRelativeError = norm(J-setup.Pathway.JBaseline,'fro') / ...
            max(norm(setup.Pathway.JBaseline,'fro'),eps);
    end
    if spec.isRepresentative
        eigenvaluesTrue = eig(full(J),'vector');
    end
else
    moved = fixedResult.State;
end

n = numel(baseline)/3;
hc = local_hcnorm(moved,baseline,n,context.CWeight);
relativeShift = norm(moved-baseline)/max(norm(baseline),eps);
caseTag = sprintf('%03d_%s_beta_%s',taskId,lower(char(spec.pathway)), ...
    local_number_tag(spec.beta));
caseDir = fullfile(outputRoot,caseTag);
if ~exist(caseDir,'dir'); mkdir(caseDir); end

summary = table(taskId,spec.pathway,spec.beta,spec.pathwayGain,-spec.beta, ...
    gain6,gainI,fixedResult.Converged,string(fixedResult.Termination), ...
    fixedResult.Iterations,fixedResult.Residual,relativeShift,hc, ...
    fixedResult.MinimumRateHz,fixedResult.MaximumRateHz,maxRealTrue, ...
    leadingImagTrue,maxRealFpp,leadingImagFpp,maxRealTrue-maxRealFpp, ...
    jvpRelativeError,baselineJacobianRelativeError,spec.isRepresentative, ...
    'VariableNames',{'taskId','pathway','beta','pathwayGain','matchedFppW', ...
    'gain6','gainI','fixedPointConverged','termination','iterations', ...
    'fixedPointResidual','relativeFixedPointShift','HCnormFromBaseline', ...
    'minimumRateHz','maximumRateHz','trueMaxRealLambda','trueLeadingImag', ...
    'fppMaxRealLambda','fppLeadingImag','trueMinusFppMaxReal', ...
    'analyticJvpRelativeError','baselineJacobianRelativeError','isRepresentative'});
writetable(summary,fullfile(caseDir,'case_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
metadata = struct('Spec',spec,'Gain6',gain6,'GainI',gainI, ...
    'MatchedFppW',-spec.beta,'FixedPointResult',fixedResult, ...
    'AnalyticJvpRelativeError',jvpRelativeError);
save(fullfile(caseDir,'case_result.mat'),'summary','metadata','moved', ...
    'eigenvaluesTrue','-v7.3');
fprintf(['Completed %s: converged=%d, residual=%.3e, shift=%.4g, HC=%.4g, ' ...
    'maxRe(true/FPP)=%.6g/%.6g, JVP error=%.3e.\n'],caseTag, ...
    fixedResult.Converged,fixedResult.Residual,relativeShift,hc, ...
    maxRealTrue,maxRealFpp,jvpRelativeError);
end

function [maximumReal,leadingImag] = local_leading(matrix)
options = struct('tol',1e-9,'maxit',1800,'p',80,'isreal',true,'disp',0);
try
    values = eigs(matrix,8,'largestreal',options);
catch
    values = eigs(matrix,8,'lr',options);
end
[maximumReal,index] = max(real(values));
leadingImag = imag(values(index));
end

function errorValue = local_jvp_check(J,phi,state,seed)
rng(4100+seed,'twister');
errors = nan(2,1);
stateScale = max(1,norm(state)/sqrt(numel(state)));
for index = 1:2
    direction = randn(size(state));
    direction = direction/norm(direction);
    step = 2e-6*stateScale;
    numeric = (phi(state+step*direction)-phi(state-step*direction))/(2*step);
    analytic = J*direction;
    errors(index) = norm(analytic-numeric)/max(norm(numeric),eps);
end
errorValue = max(errors);
end

function value = local_hcnorm(state,baseline,n,wC)
dS = state(1:n)-baseline(1:n);
dC = state(n+(1:n))-baseline(n+(1:n));
dI = state(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-wC)*dS+wC*dC;
value = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end

function tag = local_number_tag(value)
tag = sprintf('%+.3f',value);
tag = strrep(tag,'+','p');
tag = strrep(tag,'-','m');
tag = strrep(tag,'.','p');
end
