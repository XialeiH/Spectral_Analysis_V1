function result = compare_hc5_snapshot_finite_time_task(taskId,setupFile,trajectoryFile,outputRoot)
% Compare one saved HC=5 snapshot mode with its future tangent propagator.

arguments
    taskId (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    trajectoryFile (1,:) char
    outputRoot (1,:) char
end

requestedTimes = [3 15 75 78 150];
if taskId>numel(requestedTimes)
    error('Task ID must be between 1 and %d.',numel(requestedTimes));
end
startTimeMs = requestedTimes(taskId);

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
dimension = 3*n;
tauMs = setup.Config.TauMs;
l6Weight = 0;

trajectoryLoaded = load(trajectoryFile,'output');
savedSnapshots = trajectoryLoaded.output.AllSnapshots;
savedTimes = [savedSnapshots.TimeMs];
[timeError,snapshotIndex] = min(abs(savedTimes-startTimeMs));
if timeError>1e-9
    error('The requested %.3f ms snapshot is absent.',startTimeMs);
end
startState = savedSnapshots(snapshotIndex).State(:);
startHCnorm = savedSnapshots(snapshotIndex).HCnorm;

jSnapshot = sparse(l6ns_state_jacobian(setup,'L6',startState,l6Weight));
[snapshotOutput,snapshotSigma,snapshotInput] = svds(jSnapshot,1,'largest');
snapshotInput = local_orient_e(snapshotInput,n,context.CWeight);
if real(snapshotOutput'*snapshotInput)<0
    snapshotOutput = -snapshotOutput;
end

horizonGrid = [0.5 1 2 3 5 8 12 15 20 30 45 60 80 100];
stepMs = 0.5;
relativeTimeGrid = 0:stepMs:max(horizonGrid);
phi = @(state) l6ns_phi(state,l6Weight,context,[1 1],1,'fpp');
rhs = @(~,state) (phi(state)-state)/tauMs;
options = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);
solution = ode45(rhs,[0 max(relativeTimeGrid)],startState,options);
states = deval(solution,relativeTimeGrid);

identity = speye(dimension);
stepMatrices = cell(numel(relativeTimeGrid)-1,1);
for index = 1:numel(stepMatrices)
    jacobian = sparse(l6ns_state_jacobian( ...
        setup,'L6',states(:,index),l6Weight));
    stepMatrices{index} = identity + ...
        (stepMs/tauMs)*(jacobian-identity);
    if mod(index,40)==0 || index==numel(stepMatrices)
        fprintf('Task %d, t0 %.1f: tangent step %d/%d.\n', ...
            taskId,startTimeMs,index,numel(stepMatrices));
    end
end

finiteGain = zeros(size(horizonGrid));
finiteInput = cell(size(horizonGrid));
finiteOutput = cell(size(horizonGrid));
powerIterations = zeros(size(horizonGrid));
seed = snapshotInput;
for index = 1:numel(horizonGrid)
    stepCount = round(horizonGrid(index)/stepMs);
    [finiteInput{index},finiteOutput{index},finiteGain(index),powerIterations(index)] = ...
        local_top_product_singular(stepMatrices,stepCount,seed);
    if real(finiteInput{index}'*snapshotInput)<0
        finiteInput{index} = -finiteInput{index};
        finiteOutput{index} = -finiteOutput{index};
    end
    seed = finiteInput{index};
end

[peakGain,peakIndex] = max(finiteGain);
peakTimeMs = horizonGrid(peakIndex);
optimalInput = finiteInput{peakIndex};
optimalOutput = finiteOutput{peakIndex};
inputCosine = abs(real(snapshotInput'*optimalInput)) / ...
    max(norm(snapshotInput)*norm(optimalInput),eps);
outputCosine = abs(real(snapshotOutput'*optimalOutput)) / ...
    max(norm(snapshotOutput)*norm(optimalOutput),eps);

result = struct('TaskId',taskId,'StartTimeMs',startTimeMs, ...
    'StartHCnorm',startHCnorm,'SnapshotSigma',snapshotSigma, ...
    'SnapshotInput',snapshotInput,'SnapshotOutput',snapshotOutput, ...
    'FiniteTimeOptimalInput',optimalInput, ...
    'FiniteTimeOptimalOutput',optimalOutput, ...
    'PeakHorizonMs',peakTimeMs,'PeakGain',peakGain, ...
    'InputAbsoluteCosine',inputCosine, ...
    'OutputAbsoluteCosine',outputCosine, ...
    'HorizonGridMs',horizonGrid,'FiniteTimeGain',finiteGain, ...
    'PowerIterations',powerIterations,'TangentStepMs',stepMs, ...
    'MapSize',mapSize,'CWeight',context.CWeight);

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
stem = sprintf('hc5_snapshot_%02d_t%05.1fms',taskId,startTimeMs);
stem = strrep(stem,'.','p');
save(fullfile(outputRoot,[stem '.mat']),'result','-v7.3');
summary = table(taskId,startTimeMs,startHCnorm,snapshotSigma,peakTimeMs, ...
    peakGain,inputCosine,outputCosine,stepMs, ...
    'VariableNames',{'taskId','startTimeMs','startHCnorm', ...
    'snapshotDPhiSigma1','peakHorizonMs','finiteTimePeakGain', ...
    'inputAbsoluteCosine','outputAbsoluteCosine','tangentStepMs'});
writetable(summary,fullfile(outputRoot,[stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
fprintf(['Task %d complete: t0 %.1f ms, HC %.6f, Tpeak %.1f ms, ' ...
    'gain %.6g, input cosine %.6f, output cosine %.6f.\n'], ...
    taskId,startTimeMs,startHCnorm,peakTimeMs,peakGain,inputCosine,outputCosine);
end

function [right,left,sigma,iteration] = ...
        local_top_product_singular(stepMatrices,stepCount,seed)
right = real(seed(:));
right = right/max(norm(right),eps);
maximumIterations = 80;
tolerance = 1e-9;
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
