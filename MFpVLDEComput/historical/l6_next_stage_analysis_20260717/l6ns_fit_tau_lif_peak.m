function result = l6ns_fit_tau_lif_peak(cfg, context)
% Fit physical tau by matching the current ODE E peak to the LIF PSTH peak.

outputDir = fullfile(cfg.OutputRoot, 'tau_fit_lif_peak');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

w = 0;
p = context.RelaxationP;
nState = numel(context.FixedPoint);
nPopulation = nState / 3;
assert(mod(nState, 3) == 0, 'Expected equal S/C/I population sizes.');
assert(numel(context.InitialState) == nState, ...
    'Initial state and fixed point must have the same size.');

fixed = context.FixedPoint(:);
rawFixedResponse = l6ns_phi(fixed, w, context);
context.FixedPointCorrection = fixed - rawFixedResponse;
phi = @(state) l6ns_phi(state, w, context);
fixedResidual = norm(phi(fixed)-fixed) / max(norm(fixed), eps);

initialState = context.InitialState(:);
maxIterationStep = 30;
iterationStep = (0:maxIterationStep)';
iterationE = nan(size(iterationStep));
iterationState = initialState;
for stepIndex = 1:numel(iterationStep)
    iterationE(stepIndex) = local_e_population(iterationState, ...
        nPopulation, context.CWeight);
    if stepIndex < numel(iterationStep)
        response = phi(iterationState);
        iterationState = iterationState + p * (response-iterationState);
    end
end
[iterationPeakE, iterationPeakIndex] = max(iterationE);
iterationPeakStep = iterationStep(iterationPeakIndex);
assert(iterationPeakIndex > 1 && iterationPeakIndex < numel(iterationE), ...
    'Iteration E peak is not interior to the sampled interval.');

% For tau=1, intrinsic time s and physical time satisfy t=tau*s.
intrinsicStep = 0.01;
intrinsicTimeRequested = (0:intrinsicStep:20)';
odeOptions = odeset('RelTol',2e-5,'AbsTol',1e-7,'MaxStep',0.05);
[intrinsicTime, odeState] = ode45(@(~,state) phi(state)-state, ...
    intrinsicTimeRequested, initialState, odeOptions);
odeE = nan(size(intrinsicTime));
for timeIndex = 1:numel(intrinsicTime)
    odeE(timeIndex) = local_e_population(odeState(timeIndex,:)', ...
        nPopulation, context.CWeight);
end
[~, odePeakIndex] = max(odeE);
assert(odePeakIndex > 1 && odePeakIndex < numel(odeE), ...
    'ODE E peak is not interior to the intrinsic-time interval.');
[odeIntrinsicPeakTime, odePeakE] = local_parabolic_peak( ...
    intrinsicTime, odeE, odePeakIndex);

lifPeakTimeMs = cfg.LifPeakTimeMs;
tauMs = lifPeakTimeMs / odeIntrinsicPeakTime;
iterationEpochMs = p * tauMs;
physicalOdeTimeMs = tauMs * intrinsicTime;
iterationTimeMs = iterationStep * iterationEpochMs;
iterationPeakTimeMs = iterationPeakStep * iterationEpochMs;
peakTimeMismatchMs = iterationPeakTimeMs - lifPeakTimeMs;

iterationTable = table(iterationStep, iterationTimeMs, iterationE, ...
    'VariableNames', {'iterationStep','timeMs','EPopulationRate'});
odeTable = table(intrinsicTime, physicalOdeTimeMs, odeE, ...
    'VariableNames', {'intrinsicTimeAtTau1','timeMs','EPopulationRate'});
writetable(iterationTable, fullfile(outputDir, 'iteration_E_trajectory.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(odeTable, fullfile(outputDir, 'ode_E_trajectory.tsv'), ...
    'FileType','text','Delimiter','\t');

summary = table(cfg.Angle, cfg.Contrast, w, p, cfg.LifPsthBinMs, ...
    cfg.LifPeakWindow, lifPeakTimeMs, odeIntrinsicPeakTime, tauMs, ...
    iterationEpochMs, iterationPeakStep, iterationPeakTimeMs, ...
    peakTimeMismatchMs, odePeakE, iterationPeakE, fixedResidual, ...
    string(cfg.LifReferenceFile), ...
    'VariableNames', {'angle','contrast','w','relaxationP','lifPsthBinMs', ...
    'lifPeakWindow','lifPeakTimeMs','odeIntrinsicPeakTimeAtTau1', ...
    'tauMs','iterationEpochMs','iterationPeakStep','iterationPeakTimeMs', ...
    'iterationMinusLifPeakMs','odePeakE','iterationPeakE', ...
    'fixedPointResidual','lifReferenceFile'});
writetable(summary, fullfile(outputDir, 'tau_lif_fit_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[100 100 980 560]);
hold on;
stairs(iterationTimeMs, iterationE, '-', 'Color',[0 0.4470 0.7410], ...
    'LineWidth',2, 'DisplayName','relaxed iteration');
plot(physicalOdeTimeMs, odeE, '-', 'Color',[0.8500 0.3250 0.0980], ...
    'LineWidth',2, 'DisplayName',sprintf('ODE, \\tau=%.3f ms',tauMs));
xline(lifPeakTimeMs, 'k--', sprintf('LIF peak %.1f ms',lifPeakTimeMs), ...
    'LabelVerticalAlignment','bottom','DisplayName','LIF E peak');
plot(lifPeakTimeMs, odePeakE, 'o', ...
    'MarkerFaceColor',[0.8500 0.3250 0.0980], 'MarkerEdgeColor','none', ...
    'HandleVisibility','off');
plot(iterationPeakTimeMs, iterationPeakE, 'o', ...
    'MarkerFaceColor',[0 0.4470 0.7410], 'MarkerEdgeColor','none', ...
    'HandleVisibility','off');
grid on; box on;
xlabel('time (ms)');
ylabel('mean E firing rate');
title(sprintf(['LIF-calibrated ODE and iteration, angle %.1f^\\circ, ', ...
    'contrast %g, w=0'],cfg.Angle,cfg.Contrast));
legend('Location','best');
l6ns_save_figure(fig, outputDir, 'lif_calibrated_ode_iteration_E_trajectories');
close(fig);

result = struct('Summary',summary,'Iteration',iterationTable,'ODE',odeTable, ...
    'InitialState',initialState,'FixedPoint',fixed);
save(fullfile(outputDir,'tau_lif_fit_result.mat'),'result','-v7.3');
fprintf(['LIF peak %.6g ms; ODE intrinsic peak %.12g; fitted tau ', ...
    '%.12g ms; iteration epoch %.12g ms; iteration peak %.12g ms.\n'], ...
    lifPeakTimeMs,odeIntrinsicPeakTime,tauMs,iterationEpochMs, ...
    iterationPeakTimeMs);
end

function value = local_e_population(state, nPopulation, cWeight)
s = state(1:nPopulation);
c = state(nPopulation+(1:nPopulation));
value = mean((1-cWeight)*s + cWeight*c);
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
