function result = l6ns_match_tau_iteration_ode(cfg, data, context)
% Match tau in tau*df/dt=Phi(f)-f to the relaxed-iteration peak time.

outputDir = fullfile(cfg.OutputRoot, 'tau_match_iteration_ode');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

w = 0;
p = context.RelaxationP;
fixed = context.FixedPoint(:);
nState = numel(fixed);
nPopulation = nState / 3;
assert(mod(nState, 3) == 0, 'Expected equal S/C/I population sizes.');

rawFixedResponse = l6ns_phi(fixed, w, context);
context.FixedPointCorrection = fixed - rawFixedResponse;
phi = @(state) l6ns_phi(state, w, context);
fixedResidual = norm(phi(fixed) - fixed) / max(norm(fixed), eps);

populationScale = [norm(fixed(1:nPopulation))/sqrt(nPopulation), ...
    norm(fixed(nPopulation+(1:nPopulation)))/sqrt(nPopulation), ...
    norm(fixed(2*nPopulation+(1:nPopulation)))/sqrt(nPopulation)];
populationScale = max(populationScale, 1e-8);
weightVector = [ones(nPopulation,1)/populationScale(1); ...
    ones(nPopulation,1)/populationScale(2); ...
    ones(nPopulation,1)/populationScale(3)];
weightMatrix = spdiags(weightVector, 0, nState, nState);
inverseWeightMatrix = spdiags(1./weightVector, 0, nState, nState);

% J is D Phi.  The actual iteration uses M=(1-p)I+pJ.
jPhi = data.A + data.B;
jWeighted = weightMatrix * jPhi * inverseWeightMatrix;
mWeighted = (1-p) * speye(nState) + p * jWeighted;

% Use the actual relaxed iteration's one-step optimal input.  This is the
% same perturbation for both nonlinear trajectories.
rng(cfg.RandomSeed, 'twister');
inputWeighted = randn(nState, 1);
inputWeighted = inputWeighted / norm(inputWeighted);
for powerIteration = 1:30
    outputWeighted = mWeighted * inputWeighted;
    oneStepGain = norm(outputWeighted);
    outputWeighted = outputWeighted / max(oneStepGain, eps);
    inputNew = mWeighted' * outputWeighted;
    inputNew = inputNew / max(norm(inputNew), eps);
    if abs(inputNew' * inputWeighted) > 1-1e-11
        inputWeighted = inputNew;
        break
    end
    inputWeighted = inputNew;
end
oneStepGain = norm(mWeighted * inputWeighted);

requestedInitialRms = 1e-3;
physicalDirection = inverseWeightMatrix * inputWeighted;
amplitude = requestedInitialRms * sqrt(nState);
negativeDirection = physicalDirection < 0;
if any(negativeDirection)
    positiveStateLimit = 0.25 * min(fixed(negativeDirection) ./ ...
        (-physicalDirection(negativeDirection)));
    amplitude = min(amplitude, positiveStateLimit);
end
assert(isfinite(amplitude) && amplitude > 0, 'Invalid perturbation amplitude.');
initialState = fixed + amplitude * physicalDirection;
initialWeightedRms = local_weighted_rms(initialState-fixed, weightVector);

maxIterationStep = 30;
iterationStep = (0:maxIterationStep)';
iterationGain = nan(size(iterationStep));
iterationState = initialState;
for stepIndex = 1:numel(iterationStep)
    iterationGain(stepIndex) = local_weighted_rms(iterationState-fixed, ...
        weightVector) / initialWeightedRms;
    if stepIndex < numel(iterationStep)
        response = phi(iterationState);
        iterationState = iterationState + p * (response-iterationState);
    end
end
[iterationPeakGain, iterationPeakIndex] = max(iterationGain);
iterationPeakStep = iterationStep(iterationPeakIndex);
assert(iterationPeakStep > 0 && iterationPeakStep < maxIterationStep, ...
    'Iteration peak is not interior to the sampled interval.');

% For an autonomous ODE, the tau=1 solution at intrinsic time s becomes
% the tau solution at physical time t=tau*s.  Tune tau after one solve.
intrinsicStep = 0.025;
intrinsicTime = (0:intrinsicStep:20)';
odeOptions = odeset('RelTol',2e-5,'AbsTol',1e-7,'MaxStep',0.05);
[intrinsicTime, odeState] = ode45(@(~,state) phi(state)-state, ...
    intrinsicTime, initialState, odeOptions);
odeGain = nan(size(intrinsicTime));
for timeIndex = 1:numel(intrinsicTime)
    odeGain(timeIndex) = local_weighted_rms(odeState(timeIndex,:)'-fixed, ...
        weightVector) / initialWeightedRms;
end
[odePeakGainSampled, odePeakIndex] = max(odeGain);
assert(odePeakIndex > 1 && odePeakIndex < numel(odeGain), ...
    'ODE peak is not interior to the intrinsic-time interval.');
[odeIntrinsicPeakTime, odePeakGain] = local_parabolic_peak( ...
    intrinsicTime, odeGain, odePeakIndex);
tauMatched = iterationPeakStep / odeIntrinsicPeakTime;
physicalOdeTime = tauMatched * intrinsicTime;
odePhysicalPeakTime = tauMatched * odeIntrinsicPeakTime;

iterationTable = table(iterationStep, iterationGain, ...
    'VariableNames', {'timeInIterationSteps','normalizedPerturbation'});
odeTable = table(intrinsicTime, physicalOdeTime, odeGain, ...
    'VariableNames', {'intrinsicTimeAtTau1','timeInIterationSteps', ...
    'normalizedPerturbation'});
writetable(iterationTable, fullfile(outputDir, 'iteration_trajectory.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(odeTable, fullfile(outputDir, 'ode_trajectory.tsv'), ...
    'FileType','text','Delimiter','\t');

summary = table(cfg.Angle, cfg.Contrast, w, p, 1/p, tauMatched, ...
    iterationPeakStep, odeIntrinsicPeakTime, odePhysicalPeakTime, ...
    iterationPeakGain, odePeakGain, odePeakGainSampled, oneStepGain, ...
    initialWeightedRms, fixedResidual, powerIteration, ...
    'VariableNames', {'angle','contrast','w','relaxationP','eulerTau', ...
    'matchedTau','iterationPeakStep','odeIntrinsicPeakTimeAtTau1', ...
    'odePeakTimeAfterScaling','iterationPeakGain','odePeakGain', ...
    'odeSampledPeakGain','oneStepLinearGain','initialWeightedRms', ...
    'fixedPointResidual','singularPowerIterations'});
writetable(summary, fullfile(outputDir, 'tau_match_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[100 100 980 560]);
hold on;
stairs(iterationStep, iterationGain, '-', 'Color',[0 0.4470 0.7410], ...
    'LineWidth',2.0, 'DisplayName','relaxed iteration');
validOde = physicalOdeTime <= maxIterationStep;
plot(physicalOdeTime(validOde), odeGain(validOde), '-', ...
    'Color',[0.8500 0.3250 0.0980], 'LineWidth',2.0, ...
    'DisplayName',sprintf('ODE, \\tau=%.4f',tauMatched));
plot(iterationPeakStep, iterationPeakGain, 'o', ...
    'MarkerFaceColor',[0 0.4470 0.7410], 'MarkerEdgeColor','none', ...
    'HandleVisibility','off');
plot(odePhysicalPeakTime, odePeakGain, 'o', ...
    'MarkerFaceColor',[0.8500 0.3250 0.0980], 'MarkerEdgeColor','none', ...
    'HandleVisibility','off');
xline(iterationPeakStep, 'k--', sprintf('matched peak: %.2f',iterationPeakStep), ...
    'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
grid on; box on;
xlabel('time (iteration-step units)');
ylabel('normalized fixed-point-scaled perturbation');
title(sprintf('Nonlinear ODE and iteration trajectory, angle %.1f^\\circ, contrast %g, w=0', ...
    cfg.Angle,cfg.Contrast));
legend('Location','best');
l6ns_save_figure(fig, outputDir, 'matched_tau_ode_iteration_trajectories');
close(fig);

result = struct('Summary',summary,'Iteration',iterationTable,'ODE',odeTable, ...
    'PopulationScale',populationScale,'InitialInputWeighted',inputWeighted, ...
    'InitialState',initialState);
save(fullfile(outputDir,'tau_match_result.mat'),'result','-v7.3');
fprintf(['Matched tau %.12g: iteration peak step %.6g, ODE intrinsic peak ', ...
    '%.12g, scaled ODE peak %.6g.\n'], tauMatched,iterationPeakStep, ...
    odeIntrinsicPeakTime,odePhysicalPeakTime);
end

function value = local_weighted_rms(displacement, weightVector)
value = norm(weightVector .* displacement) / sqrt(numel(displacement));
end

function [peakTime, peakValue] = local_parabolic_peak(time, value, peakIndex)
indices = peakIndex + (-1:1);
x = time(indices);
y = value(indices);
coefficient = polyfit(x, y, 2);
if coefficient(1) >= 0
    peakTime = time(peakIndex);
    peakValue = value(peakIndex);
    return
end
peakTime = -coefficient(2)/(2*coefficient(1));
peakTime = min(max(peakTime,x(1)),x(end));
peakValue = polyval(coefficient,peakTime);
end
