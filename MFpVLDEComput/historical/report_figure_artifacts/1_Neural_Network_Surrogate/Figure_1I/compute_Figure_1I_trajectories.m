function compute_Figure_1I_trajectories()
% Compute baseline and two Figure 5E compensation return trajectories.

artifactRoot = fileparts(mfilename('fullpath'));
projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project';
mainRoot = fullfile(projectRoot,'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model','Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main');
addpath(fullfile(mainRoot,'Utils'),'-begin');
addpath(artifactRoot,'-begin');

sourceFile = fullfile(projectRoot,'tmp','pdfs', ...
    'perturbation7_return.qfJHU8', ...
    '7.5_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC1_data.mat');
dataFile = fullfile(mainRoot,'Figures', ...
    'spectral_analysis_eigenvalue_eigenvectors','L6 and Inhibition', ...
    'L6_and_Inhibition_offset','Figure6_0_FourContrast_FiveMetric', ...
    'figure6_0_contrast100_dataset.mat');

source = load(sourceFile,'output');
dataset = load(dataFile,'setup','betaGrid','fixedPointStates');
baseline = double(dataset.setup.BaselineState(:));
initialDisplacement = double(source.output.InitialState(:)) - baseline;
tauMs = double(source.output.TauMs);
timesMs = 0:180;

targetBeta6 = [0 0.075 0.30];
targetBetaI = 0.561654*targetBeta6;
betaGrid = double(dataset.betaGrid(:).');
actualBeta6 = zeros(1,3);
actualBetaI = zeros(1,3);
fixedPoints = zeros(numel(baseline),3);
fixedPoints(:,1) = baseline;
for conditionIndex = 2:3
    [~,beta6Index] = min(abs(betaGrid-targetBeta6(conditionIndex)));
    [~,betaIIndex] = min(abs(betaGrid-targetBetaI(conditionIndex)));
    pointIndex = betaIIndex + (beta6Index-1)*numel(betaGrid);
    actualBeta6(conditionIndex) = betaGrid(beta6Index);
    actualBetaI(conditionIndex) = betaGrid(betaIIndex);
    fixedPoints(:,conditionIndex) = double(dataset.fixedPointStates(:,pointIndex));
end

trajectories = zeros(numel(timesMs),3);
states = fixedPoints + initialDisplacement;
for conditionIndex = 1:3
    context = dataset.setup.Context;
    beta6 = actualBeta6(conditionIndex);
    betaI = actualBetaI(conditionIndex);
    phi = @(state)l6ns_phi(state,-beta6,context,[1 1],1+betaI,'true');
    rhs = @(~,state)(phi(state)-state)/tauMs;
    options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1);
    solution = ode45(rhs,[0 timesMs(end)],states(:,conditionIndex),options);
    conditionStates = deval(solution,timesMs);
    for timeIndex = 1:numel(timesMs)
        trajectories(timeIndex,conditionIndex) = local_hc( ...
            conditionStates(:,timeIndex),fixedPoints(:,conditionIndex), ...
            context.CWeight);
    end
end

cachedMask = source.output.Times<=180;
cachedBaseline = double(source.output.HCnorm(cachedMask)).';
baselineMaximumDifference = max(abs(trajectories(:,1)-cachedBaseline));
baselineRelativeRms = rms(trajectories(:,1)-cachedBaseline) / ...
    max(rms(cachedBaseline),eps);
fprintf('Baseline trajectory validation: max abs %.6g sp/s, relative RMS %.6g.\n', ...
    baselineMaximumDifference,baselineRelativeRms);
fprintf('Conditions: baseline; beta6 %.3f betaI %.3f; beta6 %.3f betaI %.3f.\n', ...
    actualBeta6(2),actualBetaI(2),actualBeta6(3),actualBetaI(3));

save(fullfile(artifactRoot,'Figure_1I_trajectories.mat'), ...
    'timesMs','trajectories','actualBeta6','actualBetaI', ...
    'baselineMaximumDifference','baselineRelativeRms','tauMs');
end

function value = local_hc(state,fixedPoint,cWeight)
n = numel(state)/3;
delta = state(:)-fixedPoint(:);
deltaE = (1-cWeight)*delta(1:n) + cWeight*delta(n+(1:n));
deltaI = delta(2*n+(1:n));
value = sqrt(0.8*mean(deltaE.^2) + 0.2*mean(deltaI.^2));
end
