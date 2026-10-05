function run_top_singular_grid_task()
% Evaluate one pathway-pair point using leading-input-mode HCnorm.

taskId = local_task_id();
gridStep = 0.008;
gridMaximum = 0.4;
gridCount = round(gridMaximum/gridStep)+1;
totalCount = gridCount^2;
if ~isfinite(taskId) || taskId<1 || taskId>totalCount || taskId~=round(taskId)
    error('TopSingularHC:TaskId','Task id must be an integer in 1:%d.',totalCount);
end

setupFile = getenv('TSHC_SETUP_FILE');
baselineModeFile = getenv('TSHC_BASELINE_MODE_FILE');
outputRoot = getenv('TSHC_OUTPUT_ROOT');
pairType = string(getenv('TSHC_PAIR'));
fixedTolerance = str2double(getenv('TSHC_FIXED_TOLERANCE'));
if ~isfile(setupFile) || ~isfile(baselineModeFile) || isempty(outputRoot)
    error('TopSingularHC:Environment','Required setup, mode, or output path is missing.');
end
if ~isfinite(fixedTolerance); fixedTolerance=1e-4; end
if ~any(pairType==["l6_inhibition","l6_excitation","excitation_inhibition"])
    error('TopSingularHC:PairType','Unsupported pair type: %s.',pairType);
end

loaded = load(setupFile,'setup');
context = loaded.setup.Context;
baseline = context.FixedPoint(:);
modeLoaded = load(baselineModeFile,'rightMode','singularValues');
baselineRightMode = modeLoaded.rightMode(:);
baselineSigma1 = modeLoaded.singularValues(1);
if numel(baselineRightMode)~=numel(baseline) || abs(norm(baselineRightMode)-1)>1e-8
    error('TopSingularHC:BaselineMode','Baseline mode has the wrong size or norm.');
end
if abs(context.CWeight-0.3077)>1e-12
    error('TopSingularHC:HCWeight','Canonical C weight is not 0.3077.');
end

betaXIndex = floor((taskId-1)/gridCount)+1;
betaYIndex = mod(taskId-1,gridCount)+1;
betaX = (betaXIndex-1)*gridStep;
betaY = (betaYIndex-1)*gridStep;
beta6 = 0;
betaE = 0;
betaI = 0;
switch pairType
    case "l6_inhibition"
        beta6 = betaX;
        betaI = betaY;
    case "l6_excitation"
        betaY = -gridMaximum+(betaYIndex-1)*gridStep;
        beta6 = betaX;
        betaE = betaY;
    case "excitation_inhibition"
        betaE = betaX;
        betaI = betaY;
end
gain6 = 1+beta6;
gainE = 1+betaE;
gainI = 1+betaI;

if pairType=="l6_inhibition"
    phi = @(state)l6ns_phi(state,1-gain6,context,[1 1],gainI,'true');
else
    phi = @(state)l6ns_phi_three_pathway(state,gain6,gainE,gainI,context);
end
fixed = real_tuning_fixed_point_tolerance(phi,baseline, ...
    context.RelaxationP,fixedTolerance);

sigma1 = NaN;
sigma2 = NaN;
gapRelative = NaN;
modeHCnorm = NaN;
modeOverlap = NaN;
singularResidual = NaN;
svdsSeconds = NaN;
if fixed.Converged
    if pairType=="l6_inhibition"
        J = real_tuning_true_jacobian(fixed.State,context,gain6,gainI);
    else
        J = real_tuning_true_jacobian_three_pathway( ...
            fixed.State,context,gain6,gainE,gainI);
    end
    options = struct('tol',1e-7,'maxit',1500,'p',60,'disp',0, ...
        'v0',baselineRightMode);
    timer = tic;
    [leftVectors,singularValues,rightVectors] = svds(J,2,'largest',options);
    svdsSeconds = toc(timer);
    singularValues = diag(singularValues);
    [singularValues,order] = sort(real(singularValues),'descend');
    leftVectors = leftVectors(:,order);
    rightVectors = rightVectors(:,order);
    sigma1 = singularValues(1);
    sigma2 = singularValues(2);
    gapRelative = (sigma1-sigma2)/max(sigma1,eps);
    [rightMode,leftMode] = local_real_singular_pair( ...
        rightVectors(:,1),leftVectors(:,1));
    modeOverlap = dot(baselineRightMode,rightMode);
    if modeOverlap<0
        rightMode = -rightMode;
        leftMode = -leftMode;
        modeOverlap = -modeOverlap;
    end
    singularResidual = max( ...
        norm(J*rightMode-sigma1*leftMode), ...
        norm(J'*leftMode-sigma1*rightMode))/max(sigma1,eps);
    n = numel(rightMode)/3;
    modeHCnorm = HC_norm_diff( ...
        rightMode(1:n),rightMode(n+(1:n)),rightMode(2*n+(1:n)), ...
        baselineRightMode(1:n),baselineRightMode(n+(1:n)), ...
        baselineRightMode(2*n+(1:n)),0.3077,0.8,0.2);
end

pointsRoot = fullfile(outputRoot,'points');
if ~exist(pointsRoot,'dir'); mkdir(pointsRoot); end
outputFile = fullfile(pointsRoot,sprintf('point_%05d.tsv',taskId));
temporaryFile = [outputFile '.tmp'];
fileId = fopen(temporaryFile,'w');
if fileId<0; error('TopSingularHC:Output','Cannot open %s.',temporaryFile); end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,['%d\t%d\t%d\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t' ...
    '%d\t%.15g\t%d\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t' ...
    '%.15g\t%.15g\t%.15g\n'],taskId,betaXIndex,betaYIndex,betaX,betaY, ...
    beta6,betaE,betaI,fixed.Converged,fixed.Residual,fixed.Iterations, ...
    modeHCnorm,sigma1,sigma2,gapRelative,modeOverlap,singularResidual, ...
    svdsSeconds,fixed.MinimumRateHz,fixed.MaximumRateHz);
clear cleanup
movefile(temporaryFile,outputFile,'f');
fprintf(['Top singular HC %s %d/%d: beta6=%.3f betaE=%.3f betaI=%.3f ' ...
    'fixed=%d residual=%.3g sigma1=%.6g baselineSigma1=%.6g HC=%.6g.\n'], ...
    pairType,taskId,totalCount,beta6,betaE,betaI,fixed.Converged, ...
    fixed.Residual,sigma1,baselineSigma1,modeHCnorm);
end

function taskId = local_task_id()
localTaskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isfinite(localTaskId)
    offset = str2double(getenv('TSHC_TASK_OFFSET'));
    if ~isfinite(offset); offset=0; end
    taskId = localTaskId+offset;
else
    taskId = str2double(getenv('TSHC_TASK_ID'));
end
end

function [rightVector,leftVector] = local_real_singular_pair(rightVector,leftVector)
[~,index] = max(abs(rightVector));
phase = exp(-1i*angle(rightVector(index)));
rightVector = real(rightVector*phase);
leftVector = real(leftVector*phase);
rightVector = rightVector/max(norm(rightVector),eps);
leftVector = leftVector/max(norm(leftVector),eps);
end
