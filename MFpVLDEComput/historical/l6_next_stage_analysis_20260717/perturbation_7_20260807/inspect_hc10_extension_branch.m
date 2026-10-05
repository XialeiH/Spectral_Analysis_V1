function report = inspect_hc10_extension_branch(setupFile,longDataFile,correctedRuntime,outputRoot)
% Determine whether the HC=10 equilibrium occupies the high-L4E artifact branch.

arguments
    setupFile (1,:) char
    longDataFile (1,:) char
    correctedRuntime (1,:) char
    outputRoot (1,:) char
end

setupLoaded = load(setupFile,'setup');
setup = setupLoaded.setup;
context = setup.Context;
longLoaded = load(longDataFile,'result');
highState = longLoaded.result.States(:,end);
originalState = context.FixedPoint(:);
n = numel(originalState)/3;

originalInputs = local_inputs(originalState,context);
highInputs = local_inputs(highState,context);
trainingEdge = 47500;
guardFull = 55000;

inputName = ["L4E_S";"L4E_C";"L4E_I"];
originalStats = local_stats(originalInputs.L4E,trainingEdge,guardFull);
highStats = local_stats(highInputs.L4E,trainingEdge,guardFull);
inputSummary = table(inputName, ...
    originalStats.minimum,originalStats.median,originalStats.maximum, ...
    originalStats.fractionAboveEdge,originalStats.fractionAboveFull, ...
    highStats.minimum,highStats.median,highStats.maximum, ...
    highStats.fractionAboveEdge,highStats.fractionAboveFull, ...
    'VariableNames',{'target','originalMin','originalMedian','originalMax', ...
    'originalFractionAbove47500','originalFractionAbove55000', ...
    'highMin','highMedian','highMax','highFractionAbove47500', ...
    'highFractionAbove55000'});

phiUnguarded = @(state) l6ns_phi(state,0,context,[1 1],1,'fpp');
unguardedResidualHigh = local_hcnorm(phiUnguarded(highState),highState,n,context.CWeight);

addpath(correctedRuntime,'-begin');
clear LocalResponse_6D_MLP_prefAngle
resolvedLocalResponse = which('LocalResponse_6D_MLP_prefAngle');
phiGuarded = @(state) l6ns_phi(state,0,context,[1 1],1,'fpp');
guardedResidualOriginal = local_hcnorm(phiGuarded(originalState),originalState,n,context.CWeight);
guardedMappedHigh = phiGuarded(highState);
guardedResidualHigh = local_hcnorm(guardedMappedHigh,highState,n,context.CWeight);
guardedOneStepHCFromOriginal = local_hcnorm(guardedMappedHigh,originalState,n,context.CWeight);

guardedState = highState;
iterationHCFromOriginal = zeros(100,1);
iterationResidualHC = zeros(100,1);
for index = 1:100
    mapped = phiGuarded(guardedState);
    iterationResidualHC(index) = local_hcnorm(mapped,guardedState,n,context.CWeight);
    guardedState = mapped;
    iterationHCFromOriginal(index) = local_hcnorm(guardedState,originalState,n,context.CWeight);
end
guardedFinalResidual = local_hcnorm(phiGuarded(guardedState),guardedState,n,context.CWeight);

guardComparison = table(unguardedResidualHigh,guardedResidualOriginal, ...
    guardedResidualHigh,guardedOneStepHCFromOriginal, ...
    iterationHCFromOriginal(end),guardedFinalResidual, ...
    string(resolvedLocalResponse), ...
    'VariableNames',{'unguardedResidualAtHighState', ...
    'guardedResidualAtOriginalState','guardedResidualAtHighState', ...
    'guardedOneStepHCFromOriginal','guardedIteration100HCFromOriginal', ...
    'guardedIteration100ResidualHC','resolvedGuardedLocalResponse'});

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
writetable(inputSummary,fullfile(outputRoot,'hc10_extension_input_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(guardComparison,fullfile(outputRoot,'hc10_guard_comparison.tsv'), ...
    'FileType','text','Delimiter','\t');
report = struct('InputSummary',inputSummary,'GuardComparison',guardComparison, ...
    'OriginalInputs',originalInputs,'HighInputs',highInputs, ...
    'GuardedIterationHCFromOriginal',iterationHCFromOriginal, ...
    'GuardedIterationResidualHC',iterationResidualHC, ...
    'GuardedFinalState',guardedState,'TrainingEdge',trainingEdge, ...
    'GuardFull',guardFull);
save(fullfile(outputRoot,'hc10_extension_branch_report.mat'),'report','-v7.3');

fprintf('High-state L4E fractions above 47500: S %.3f, C %.3f, I %.3f.\n', ...
    highStats.fractionAboveEdge);
fprintf(['Unguarded residual %.9g; guarded high-state residual %.9g; ' ...
    'after 100 guarded iterations HC(original) %.9g, residual %.9g.\n'], ...
    unguardedResidualHigh,guardedResidualHigh, ...
    iterationHCFromOriginal(end),guardedFinalResidual);
end

function inputs = local_inputs(state,context)
n = numel(state)/3;
s = state(1:n);
c = state(n+(1:n));
i = state(2*n+(1:n));
if context.Isaturation
    eRaw = (1-context.CWeight)*s+context.CWeight*c;
    eBase = L6Convert(eRaw,context.EKpUse);
    adjustment = eBase./eRaw;
    adjustment(~isfinite(adjustment)) = 1;
    sUse = s.*adjustment;
    cUse = c.*adjustment;
    iUse = InhMulp(i,context.IKpUse);
else
    sUse = s;
    cUse = c;
    iUse = i;
end
inputs.L4E = {(context.C_SS*sUse+context.C_SC*cUse)/context.L4SEp; ...
    (context.C_CS*sUse+context.C_CC*cUse)/context.L4CEp; ...
    (context.C_IS*sUse+context.C_IC*cUse)/context.L4IEp};
inputs.L4I = {(context.C_SI*iUse)/context.L4SIp; ...
    (context.C_CI*iUse)/context.L4CIp; ...
    (context.C_II*iUse)/context.L4IIp};
inputs.L6 = l6ns_l6_dynamic(state,context,[1 1]);
end

function stats = local_stats(values,edge,full)
count = numel(values);
stats.minimum = zeros(count,1);
stats.median = zeros(count,1);
stats.maximum = zeros(count,1);
stats.fractionAboveEdge = zeros(count,1);
stats.fractionAboveFull = zeros(count,1);
for index = 1:count
    current = values{index}(:);
    stats.minimum(index) = min(current);
    stats.median(index) = median(current);
    stats.maximum(index) = max(current);
    stats.fractionAboveEdge(index) = mean(current>edge);
    stats.fractionAboveFull(index) = mean(current>full);
end
end

function value = local_hcnorm(state,baseline,n,wC)
dS = state(1:n)-baseline(1:n);
dC = state(n+(1:n))-baseline(n+(1:n));
dI = state(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-wC)*dS+wC*dC;
value = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end
