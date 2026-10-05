function run_figure6_17_response_tangent_task()
% Finite-difference one visual-condition block of T = D_beta R(0).

taskId = str2double(getenv('FIG617_TASK_ID'));
arrayId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isfinite(arrayId); taskId = arrayId; end
contrasts = [19 42 66 100];
angleCount = 7;
if ~isfinite(taskId) || taskId ~= round(taskId) || taskId < 1 || ...
        taskId > numel(contrasts)*angleCount
    error('Figure617:TaskId','Task id must be an integer in 1:28.');
end

runRoot = getenv('FIG617_RUN_ROOT');
if isempty(runRoot); error('Figure617:Environment','FIG617_RUN_ROOT is required.'); end
step = str2double(getenv('FIG617_STEP'));
if ~isfinite(step); step = 0.008; end
tolerance = str2double(getenv('FIG617_FIXED_TOLERANCE'));
if ~isfinite(tolerance); tolerance = 1e-4; end

contrastIndex = floor((taskId-1)/angleCount)+1;
angleIndex = mod(taskId-1,angleCount)+1;
contrast = contrasts(contrastIndex);
loaded = load(fullfile(runRoot,'setup',sprintf( ...
    'figure6_0_baseline_contrast%d.mat',contrast)),'setup');
setup = loaded.setup;
angleDeg = setup.CanonicalAngles(angleIndex);
baseline = double(setup.BaselineCanonicalStates(:,angleIndex));
context = setup.Context;
context.OrientationUse = angleDeg*ones(size(context.OrientationUse));

phi6 = @(state)l6ns_phi(state,-step,context,[1 1],1,'true');
phiI = @(state)l6ns_phi(state,0,context,[1 1],1+step,'true');
fixed6 = real_tuning_fixed_point_tolerance( ...
    phi6,baseline,context.RelaxationP,tolerance);
fixedI = real_tuning_fixed_point_tolerance( ...
    phiI,baseline,context.RelaxationP,tolerance);

result = struct();
result.TaskId = taskId;
result.Contrast = contrast;
result.AngleDeg = angleDeg;
result.AngleIndex = angleIndex;
result.Step = step;
result.FixedTolerance = tolerance;
result.Beta6Accepted = fixed6.Converged;
result.BetaIAccepted = fixedI.Converged;
result.Beta6Residual = fixed6.Residual;
result.BetaIResidual = fixedI.Residual;
result.Beta6Iterations = fixed6.Iterations;
result.BetaIIterations = fixedI.Iterations;
if fixed6.Converged && fixedI.Converged
    result.dResponse_dBeta6 = single((fixed6.State-baseline)/step);
    result.dResponse_dBetaI = single((fixedI.State-baseline)/step);
else
    result.dResponse_dBeta6 = single(nan(size(baseline)));
    result.dResponse_dBetaI = single(nan(size(baseline)));
end

pointRoot = fullfile(runRoot,'results','figure6_17','points');
if ~exist(pointRoot,'dir'); mkdir(pointRoot); end
outputFile = fullfile(pointRoot,sprintf('condition_%02d.mat',taskId));
temporaryFile = [outputFile '.tmp.mat'];
save(temporaryFile,'result','-v7');
movefile(temporaryFile,outputFile,'f');
fprintf(['Figure 6.17 task %d: contrast %d angle %.2f, accepted [%d %d], ' ...
    'residual [%.3g %.3g].\n'],taskId,contrast,angleDeg, ...
    fixed6.Converged,fixedI.Converged,fixed6.Residual,fixedI.Residual);
end
