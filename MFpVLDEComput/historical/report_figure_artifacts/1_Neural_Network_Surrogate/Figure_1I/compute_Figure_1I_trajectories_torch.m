function compute_Figure_1I_trajectories_torch(sourceFile,dataFile,outputFile)
% Torch computation of the three Figure 1I nonlinear ODE trajectories.

source = load(sourceFile,'output');
dataset = load(dataFile,'setup','betaGrid','fixedPointStates');
baseline = double(dataset.setup.BaselineState(:));
initialDisplacement = double(source.output.InitialState(:))-baseline;
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
    pointIndex = betaIIndex+(beta6Index-1)*numel(betaGrid);
    actualBeta6(conditionIndex) = betaGrid(beta6Index);
    actualBetaI(conditionIndex) = betaGrid(betaIIndex);
    fixedPoints(:,conditionIndex) = double(dataset.fixedPointStates(:,pointIndex));
end

trajectories = zeros(numel(timesMs),3);
for conditionIndex = 1:3
    beta6 = actualBeta6(conditionIndex);
    betaI = actualBetaI(conditionIndex);
    context = dataset.setup.Context;
    phi = @(state)l6ns_phi(state,-beta6,context,[1 1],1+betaI,'true');
    rhs = @(~,state)(phi(state)-state)/tauMs;
    initialState = fixedPoints(:,conditionIndex)+initialDisplacement;
    options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1);
    solution = ode45(rhs,[0 timesMs(end)],initialState,options);
    states = deval(solution,timesMs);
    for timeIndex = 1:numel(timesMs)
        trajectories(timeIndex,conditionIndex) = local_hc( ...
            states(:,timeIndex),fixedPoints(:,conditionIndex),context.CWeight);
    end
    fprintf('Condition %d complete: beta6 %.3f betaI %.3f peak %.6g sp/s.\n', ...
        conditionIndex,beta6,betaI,max(trajectories(:,conditionIndex)));
end

cachedMask = source.output.Times<=180;
cachedBaseline = double(source.output.HCnorm(cachedMask)).';
baselineMaximumDifference = max(abs(trajectories(:,1)-cachedBaseline));
baselineRelativeRms = rms(trajectories(:,1)-cachedBaseline)/ ...
    max(rms(cachedBaseline),eps);
fprintf('Baseline validation: max abs %.6g sp/s, relative RMS %.6g.\n', ...
    baselineMaximumDifference,baselineRelativeRms);
save(outputFile,'timesMs','trajectories','actualBeta6','actualBetaI', ...
    'baselineMaximumDifference','baselineRelativeRms','tauMs');
end

function value = local_hc(state,fixedPoint,cWeight)
n = numel(state)/3;
delta = state(:)-fixedPoint(:);
deltaE = (1-cWeight)*delta(1:n)+cWeight*delta(n+(1:n));
deltaI = delta(2*n+(1:n));
value = sqrt(0.8*mean(deltaE.^2)+0.2*mean(deltaI.^2));
end
