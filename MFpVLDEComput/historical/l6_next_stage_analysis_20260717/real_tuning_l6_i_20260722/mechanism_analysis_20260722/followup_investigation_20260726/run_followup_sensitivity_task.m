function run_followup_sensitivity_task()
% Decompose d lambda/d beta into direct-input and equilibrium-movement terms.

paths = followup_initialize();
loaded = load(paths.SetupFile,'setup');
context = loaded.setup.Context;
specs = followup_specs('sensitivity');
taskId = local_task_id();
spec = specs(taskId,:);
beta = spec.beta;
h = 1e-4;
[gain6,gainI] = local_gains(spec.pathway,beta);
state = followup_branch_state(paths,spec.pathway,beta,"low",context);
stateMinus = followup_branch_state(paths,spec.pathway,beta-h,"low",context);
statePlus = followup_branch_state(paths,spec.pathway,beta+h,"low",context);

J = real_tuning_true_jacobian(state,context,gain6,gainI);
[lambda,right,left,eigenResidual] = followup_leading_pair(J);
[gain6Minus,gainIMinus] = local_gains(spec.pathway,beta-h);
[gain6Plus,gainIPlus] = local_gains(spec.pathway,beta+h);
JTotalMinus = real_tuning_true_jacobian( ...
    stateMinus,context,gain6Minus,gainIMinus);
JTotalPlus = real_tuning_true_jacobian( ...
    statePlus,context,gain6Plus,gainIPlus);
JDirectMinus = real_tuning_true_jacobian( ...
    state,context,gain6Minus,gainIMinus);
JDirectPlus = real_tuning_true_jacobian( ...
    state,context,gain6Plus,gainIPlus);
dJTotal = (JTotalPlus-JTotalMinus)/(2*h);
dJDirect = (JDirectPlus-JDirectMinus)/(2*h);
dJMovement = dJTotal-dJDirect;
normalizer = left'*right;
lambdaPrimeDirect = real(left'*(dJDirect*right)/normalizer);
lambdaPrimeMovement = real(left'*(dJMovement*right)/normalizer);
lambdaPrimeTotal = real(left'*(dJTotal*right)/normalizer);

lambdaMinus = local_matched_lambda(JTotalMinus,left);
lambdaPlus = local_matched_lambda(JTotalPlus,left);
lambdaPrimeFiniteDifference = real(lambdaPlus-lambdaMinus)/(2*h);

phiMinusAtState = mechanism_phi_variant( ...
    state,context,gain6Minus,gainIMinus,'extended');
phiPlusAtState = mechanism_phi_variant( ...
    state,context,gain6Plus,gainIPlus,'extended');
b = (phiPlusAtState-phiMinusAtState)/(2*h);
vFiniteDifference = (statePlus-stateMinus)/(2*h);
A = speye(numel(state))-J;
vResolvent = A\b;
vRelativeError = norm(vResolvent-vFiniteDifference)/ ...
    max(norm(vFiniteDifference),eps);
directionalResolventGain = norm(vResolvent)/max(norm(b),eps);
modeForcingProjection = left'*b/normalizer;
singleModeDisplacement = modeForcingProjection/(1-lambda);
fullResolventTwoNorm = local_resolvent_two_norm(A);

[~,~,details] = mechanism_jacobian_components( ...
    state,context,gain6,gainI);
slopeRows = local_slope_rows(taskId,spec.pathway,beta,details);
row = table(taskId,spec.pathway,beta,h,real(lambda),imag(lambda), ...
    eigenResidual,lambdaPrimeDirect,lambdaPrimeMovement,lambdaPrimeTotal, ...
    lambdaPrimeFiniteDifference, ...
    abs(lambdaPrimeTotal-lambdaPrimeFiniteDifference)/ ...
        max(abs(lambdaPrimeFiniteDifference),eps), ...
    norm(b),norm(vFiniteDifference),norm(vResolvent),vRelativeError, ...
    directionalResolventGain,fullResolventTwoNorm, ...
    real(modeForcingProjection),imag(modeForcingProjection), ...
    real(singleModeDisplacement),imag(singleModeDisplacement), ...
    'VariableNames',{'taskId','pathway','beta','finiteDifferenceStep', ...
    'leadingReal','leadingImag','eigenpairResidual', ...
    'lambdaPrimeDirect','lambdaPrimeMovement','lambdaPrimeTotal', ...
    'lambdaPrimeTrackedFiniteDifference','lambdaPrimeRelativeError', ...
    'forcingNorm','stateDerivativeFiniteDifferenceNorm','stateDerivativeResolventNorm', ...
    'stateDerivativeRelativeError','directionalResolventGain','fullResolventTwoNorm', ...
    'modeForcingProjectionReal','modeForcingProjectionImag', ...
    'singleModeDisplacementReal','singleModeDisplacementImag'});

outDir = fullfile(paths.FollowupOutput,'sensitivity_tasks');
if ~exist(outDir,'dir'); mkdir(outDir); end
writetable(row,fullfile(outDir,sprintf('sensitivity_%03d.tsv',taskId)), ...
    'FileType','text','Delimiter','\t');
writetable(slopeRows,fullfile(outDir,sprintf('slopes_%03d.tsv',taskId)), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outDir,sprintf('sensitivity_%03d.mat',taskId)), ...
    'row','slopeRows','lambda','right','left','b','vFiniteDifference', ...
    'vResolvent','-v7.3');
fprintf(['%s beta=%+.6f lambda %.9f: lambda prime direct/movement/total/FD ' ...
    '%+.6g %+.6g %+.6g %+.6g; v error %.3e; ||resolvent||2 %.4g\n'], ...
    spec.pathway,beta,real(lambda),lambdaPrimeDirect,lambdaPrimeMovement, ...
    lambdaPrimeTotal,lambdaPrimeFiniteDifference,vRelativeError, ...
    fullResolventTwoNorm);
end

function value = local_resolvent_two_norm(A)
options = struct('tol',2e-5,'maxit',1600,'p',80,'disp',0);
try
    sigmaMin = svds(A,1,'smallest',options);
    value = 1/sigmaMin;
catch exception
    warning('Followup:ResolventNorm','svds failed: %s',exception.message);
    value = NaN;
end
end

function lambda = local_matched_lambda(J,leftReference)
options = struct('tol',1e-10,'maxit',2400,'p',120,'isreal',true,'disp',0);
[vectors,values] = eigs(J,20,'largestreal',options);
values = diag(values);
scores = abs(leftReference'*vectors) ./ ...
    (norm(leftReference)*vecnorm(vectors,2,1));
[~,index] = max(scores);
lambda = values(index);
end

function rows = local_slope_rows(taskId,pathway,beta,details)
n = numel(details.L6);
targetNames = ["S";"C";"I"];
channelNames = ["L4E";"L4I";"L6"];
channels = {details.L4EGradient,details.L4IGradient,details.L6Gradient};
cells = cell(9,1);
rowIndex = 0;
for target = 1:3
    indices = (target-1)*n+(1:n);
    for channel = 1:3
        rowIndex = rowIndex+1;
        values = channels{channel}(indices);
        q = quantile(values,[0.05 0.25 0.5 0.75 0.95]);
        cells{rowIndex} = table(taskId,pathway,beta,targetNames(target), ...
            channelNames(channel),min(values),q(1),q(2),q(3),q(4),q(5), ...
            max(values),mean(abs(values)<1e-6), ...
            'VariableNames',{'taskId','pathway','beta','target','channel', ...
            'minimum','q05','q25','median','q75','q95','maximum', ...
            'fractionAbsBelow1eMinus6'});
    end
end
rows = vertcat(cells{:});
end

function taskId = local_task_id()
taskId = str2double(getenv('FOLLOWUP_TASK_ID'));
if ~isfinite(taskId); taskId=str2double(getenv('SLURM_ARRAY_TASK_ID')); end
if ~isfinite(taskId); taskId=1; end
taskId=round(taskId);
end

function [gain6,gainI] = local_gains(pathway,beta)
gain6=1; gainI=1;
if pathway=="L6"; gain6=1+beta; else; gainI=1+beta; end
end
