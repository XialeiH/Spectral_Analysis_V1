function result = compute_dense_finite_time_alignment_task(taskId,setupFile,trajectoryRoot,outputRoot,figureIndexOverride,sampleIndexOverride)
% Compute the finite-time response alignment at one 2 ms ODE sample.

arguments
    taskId (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    trajectoryRoot (1,:) char
    outputRoot (1,:) char
    figureIndexOverride (1,1) double = NaN
    sampleIndexOverride (1,1) double = NaN
end

sampleCount = 251;
figureCount = 12;
if taskId>sampleCount*figureCount
    error('Task ID %d exceeds %d.',taskId,sampleCount*figureCount);
end
figureIndex = floor((taskId-1)/sampleCount)+1;
sampleIndex = mod(taskId-1,sampleCount)+1;
explicitSample = isfinite(figureIndexOverride) && isfinite(sampleIndexOverride);
if explicitSample
    figureIndex = figureIndexOverride;
    sampleIndex = sampleIndexOverride;
end

setupLoaded = load(setupFile,'setup');
setup = setupLoaded.setup;
context = setup.Context;
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
dimension = 3*n;
tauMs = setup.Config.TauMs;
loaded = load(fullfile(trajectoryRoot, ...
    sprintf('figure_%02d_trajectory.mat',figureIndex)), ...
    'states','trajectoryTimesMs','sampleTimesMs','l6Weight', ...
    'baselineSingularOutputEMap');

sampleTimeMs = loaded.sampleTimesMs(sampleIndex);
tangentStepMs = loaded.trajectoryTimesMs(2)-loaded.trajectoryTimesMs(1);
startIndex = round(sampleTimeMs/tangentStepMs)+1;
startState = loaded.states(:,startIndex);
jSnapshot = sparse(l6ns_state_jacobian( ...
    setup,'L6',startState,loaded.l6Weight));
[~,~,snapshotInput] = svds(jSnapshot,1,'largest');
snapshotInput = local_orient_e(snapshotInput,n,context.CWeight);

shortHorizons = [0.5 1 2 3 5 8 12 15 20 30];
longHorizons = [0.5 1 2 3 5 8 12 15 20 30 45 60 80 100];
identity = speye(dimension);
stepMatrices = cell(round(max(longHorizons)/tangentStepMs),1);
shortStepCount = round(max(shortHorizons)/tangentStepMs);
for index = 1:shortStepCount
    stateIndex = startIndex+index-1;
    jacobian = sparse(l6ns_state_jacobian( ...
        setup,'L6',loaded.states(:,stateIndex),loaded.l6Weight));
    stepMatrices{index} = identity + ...
        (tangentStepMs/tauMs)*(jacobian-identity);
end
[finiteInput,finiteOutput,finiteGain] = ...
    local_scan(stepMatrices,shortHorizons,tangentStepMs,snapshotInput);

extendedTo100Ms = false;
if finiteGain(end)>=max(finiteGain(1:end-1))
    extendedTo100Ms = true;
    for index = shortStepCount+1:numel(stepMatrices)
        stateIndex = startIndex+index-1;
        jacobian = sparse(l6ns_state_jacobian( ...
            setup,'L6',loaded.states(:,stateIndex),loaded.l6Weight));
        stepMatrices{index} = identity + ...
            (tangentStepMs/tauMs)*(jacobian-identity);
    end
    [finiteInput,finiteOutput,finiteGain] = ...
        local_scan(stepMatrices,longHorizons,tangentStepMs,snapshotInput);
    horizonGrid = longHorizons;
else
    horizonGrid = shortHorizons;
end

[peakGain,peakIndex] = max(finiteGain);
optimalInput = finiteInput{peakIndex};
optimalOutput = finiteOutput{peakIndex};
if real(optimalInput'*snapshotInput)<0
    optimalOutput = -optimalOutput;
end
finiteTimeOutputEMap = local_e_map(optimalOutput,mapSize,context.CWeight);
baselineMap = loaded.baselineSingularOutputEMap;
alignment = abs(real(finiteTimeOutputEMap(:)'*baselineMap(:))) / ...
    max(norm(finiteTimeOutputEMap(:))*norm(baselineMap(:)),eps);

result = struct('TaskId',taskId,'FigureIndex',figureIndex, ...
    'SampleIndex',sampleIndex,'TimeMs',sampleTimeMs, ...
    'Alignment',alignment,'PeakHorizonMs',horizonGrid(peakIndex), ...
    'PeakGain',peakGain,'ExtendedTo100Ms',extendedTo100Ms);
if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
if explicitSample
    outputName = sprintf('figure_%02d_sample_%03d.mat',figureIndex,sampleIndex);
else
    outputName = sprintf('task_%04d.mat',taskId);
end
save(fullfile(outputRoot,outputName),'result');
fprintf('Figure %d sample %d: t=%.1f ms, cosine=%.8f.\n', ...
    figureIndex,sampleIndex,sampleTimeMs,alignment);
end

function [inputs,outputs,gains] = local_scan(stepMatrices,horizons,stepMs,seed)
inputs = cell(size(horizons));
outputs = cell(size(horizons));
gains = zeros(size(horizons));
currentSeed = seed;
for index = 1:numel(horizons)
    stepCount = round(horizons(index)/stepMs);
    [inputs{index},outputs{index},gains(index)] = ...
        local_top_product_singular(stepMatrices,stepCount,currentSeed);
    if real(inputs{index}'*seed)<0
        inputs{index} = -inputs{index};
        outputs{index} = -outputs{index};
    end
    currentSeed = inputs{index};
end
end

function [right,left,sigma] = local_top_product_singular(stepMatrices,stepCount,seed)
right = real(seed(:));
right = right/max(norm(right),eps);
maximumIterations = 40;
tolerance = 1e-8;
for iteration = 1:maximumIterations
    leftRaw = local_apply_forward(stepMatrices,stepCount,right);
    candidate = local_apply_adjoint(stepMatrices,stepCount,leftRaw);
    candidate = candidate/max(norm(candidate),eps);
    if abs(real(candidate'*right))>1-tolerance
        right = candidate;
        break
    end
    right = candidate;
end
leftRaw = local_apply_forward(stepMatrices,stepCount,right);
sigma = norm(leftRaw);
left = leftRaw/max(sigma,eps);
end

function output = local_apply_forward(stepMatrices,stepCount,input)
output = input;
for index = 1:stepCount
    output = stepMatrices{index}*output;
end
end

function output = local_apply_adjoint(stepMatrices,stepCount,input)
output = input;
for index = stepCount:-1:1
    output = stepMatrices{index}'*output;
end
end

function vector = local_orient_e(vector,n,wC)
e = (1-wC)*real(vector(1:n))+wC*real(vector(n+(1:n)));
[~,index] = max(abs(e));
if e(index)<0
    vector = -vector;
end
end

function map = local_e_map(vector,mapSize,wC)
n = prod(mapSize);
map = reshape((1-wC)*real(vector(1:n))+ ...
    wC*real(vector(n+(1:n))),mapSize);
end
