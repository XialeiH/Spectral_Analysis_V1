function result = compute_all_figure_finite_time_row_task(taskId,setupFile,sourceRoot,outputRoot)
% Compute one finite-time transient map for Figures 7.1-7.12.

arguments
    taskId (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    sourceRoot (1,:) char
    outputRoot (1,:) char
end

figureNumbers = {'7.1','7.2','7.3','7.4','7.5','7.6', ...
    '7.7','7.8','7.9','7.10','7.11','7.12'};
rowCounts = [4 4 5 5 5 6 6 6 5 5 5 5];
sourceRelativePaths = { ...
    'outputs/7.1_Perturbation_Transient_and_Return_Baseline_L6_data.mat', ...
    'outputs/7.2_Perturbation_Transient_and_Return_Enhanced_L6_1p15x_data.mat', ...
    'outputs/7.4_Perturbation_Transient_and_Return_Baseline_L6_LargePerturbation_HC0p30_data.mat', ...
    'outputs/7.5_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC1_data.mat', ...
    'outputs/7.6_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC5_data.mat', ...
    'fixed_point_row_figures/7.6_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p00_fixedpointrow_data.mat', ...
    'fixed_point_row_figures/7.7_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p01_fixedpointrow_data.mat', ...
    'fixed_point_row_figures/7.8_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p02_fixedpointrow_data.mat', ...
    'fixed_direction_boundary_figures/7.9_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p03_data.mat', ...
    'fixed_direction_boundary_figures/7.10_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p10_data.mat', ...
    'fixed_direction_boundary_figures/7.11_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p50_data.mat', ...
    'outputs/7.7_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC10_data.mat'};
if taskId>sum(rowCounts)
    error('Task ID %d exceeds the %d figure rows.',taskId,sum(rowCounts));
end
[figureIndex,rowIndex] = local_task_location(taskId,rowCounts);
figureNumber = figureNumbers{figureIndex};
sourceFile = fullfile(sourceRoot,sourceRelativePaths{figureIndex});

setupLoaded = load(setupFile,'setup');
setup = setupLoaded.setup;
context = setup.Context;
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
dimension = 3*n;
tauMs = setup.Config.TauMs;

sourceLoaded = load(sourceFile,'output');
source = sourceLoaded.output;
if isfield(source,'AllSnapshots')
    snapshots = source.AllSnapshots;
else
    snapshots = source.Snapshots;
end
if numel(snapshots)~=rowCounts(figureIndex)
    error('%s expected %d rows but contains %d.', ...
        figureNumber,rowCounts(figureIndex),numel(snapshots));
end
startState = snapshots(rowIndex).State(:);
startTimeMs = snapshots(rowIndex).TimeMs;
startHCnorm = snapshots(rowIndex).HCnorm;
l6Weight = source.L6Weight;

jSnapshot = sparse(l6ns_state_jacobian(setup,'L6',startState,l6Weight));
[~,~,snapshotInput] = svds(jSnapshot,1,'largest');
snapshotInput = local_orient_e(snapshotInput,n,context.CWeight);

stepMs = 0.5;
shortHorizons = [0.5 1 2 3 5 8 12 15 20 30];
longHorizons = [0.5 1 2 3 5 8 12 15 20 30 45 60 80 100];
relativeTimeGrid = 0:stepMs:max(longHorizons);
phi = @(state) l6ns_phi(state,l6Weight,context,[1 1],1,'fpp');
rhs = @(~,state) (phi(state)-state)/tauMs;
options = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);
solution = ode45(rhs,[0 max(relativeTimeGrid)],startState,options);
states = deval(solution,relativeTimeGrid);

identity = speye(dimension);
shortStepCount = round(max(shortHorizons)/stepMs);
stepMatrices = cell(round(max(longHorizons)/stepMs),1);
for index = 1:shortStepCount
    jacobian = sparse(l6ns_state_jacobian( ...
        setup,'L6',states(:,index),l6Weight));
    stepMatrices{index} = identity + ...
        (stepMs/tauMs)*(jacobian-identity);
end
[finiteInput,finiteOutput,finiteGain,powerIterations] = ...
    local_scan(stepMatrices,shortHorizons,stepMs,snapshotInput);

extendedTo100Ms = false;
if finiteGain(end)>=max(finiteGain(1:end-1))
    extendedTo100Ms = true;
    for index = shortStepCount+1:numel(stepMatrices)
        jacobian = sparse(l6ns_state_jacobian( ...
            setup,'L6',states(:,index),l6Weight));
        stepMatrices{index} = identity + ...
            (stepMs/tauMs)*(jacobian-identity);
    end
    [finiteInput,finiteOutput,finiteGain,powerIterations] = ...
        local_scan(stepMatrices,longHorizons,stepMs,snapshotInput);
    horizonGrid = longHorizons;
else
    horizonGrid = shortHorizons;
end

[peakGain,peakIndex] = max(finiteGain);
peakHorizonMs = horizonGrid(peakIndex);
optimalInput = finiteInput{peakIndex};
optimalOutput = finiteOutput{peakIndex};
if real(optimalInput'*snapshotInput)<0
    optimalInput = -optimalInput;
    optimalOutput = -optimalOutput;
end
inputCosine = abs(real(snapshotInput'*optimalInput)) / ...
    max(norm(snapshotInput)*norm(optimalInput),eps);
finiteTimeEMap = local_e_map(optimalInput,mapSize,context.CWeight);
finiteTimeEMap = finiteTimeEMap/max(norm(finiteTimeEMap(:)),eps);
finiteTimeOutputEMap = local_e_map(optimalOutput,mapSize,context.CWeight);

result = struct('TaskId',taskId,'FigureIndex',figureIndex, ...
    'FigureNumber',figureNumber,'RowIndex',rowIndex, ...
    'StartTimeMs',startTimeMs,'StartHCnorm',startHCnorm, ...
    'L6Weight',l6Weight,'PeakHorizonMs',peakHorizonMs, ...
    'PeakGain',peakGain,'InputAbsoluteCosine',inputCosine, ...
    'FiniteTimeOptimalInputEMap',finiteTimeEMap, ...
    'FiniteTimeOptimalOutputEMap',finiteTimeOutputEMap, ...
    'HorizonGridMs',horizonGrid,'FiniteTimeGain',finiteGain, ...
    'PowerIterations',powerIterations,'TangentStepMs',stepMs, ...
    'ExtendedTo100Ms',extendedTo100Ms);

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
stem = sprintf('figure_%02d_row_%02d',figureIndex,rowIndex);
save(fullfile(outputRoot,[stem '.mat']),'result','-v7.3');
summary = table(taskId,figureIndex,string(figureNumber),rowIndex, ...
    startTimeMs,startHCnorm,l6Weight,peakHorizonMs,peakGain,inputCosine, ...
    stepMs,extendedTo100Ms, ...
    'VariableNames',{'taskId','figureIndex','figureNumber','rowIndex', ...
    'startTimeMs','startHCnorm','l6Weight','peakHorizonMs', ...
    'finiteTimePeakGain','inputAbsoluteCosine','tangentStepMs', ...
    'extendedTo100Ms'});
writetable(summary,fullfile(outputRoot,[stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
fprintf(['%s row %d: t %.1f ms, HC %.6g, Tpeak %.1f ms, gain %.6g, ' ...
    'input cosine %.6f, extended=%d.\n'],figureNumber,rowIndex, ...
    startTimeMs,startHCnorm,peakHorizonMs,peakGain,inputCosine,extendedTo100Ms);
end

function [figureIndex,rowIndex] = local_task_location(taskId,rowCounts)
cumulative = cumsum(rowCounts);
figureIndex = find(taskId<=cumulative,1);
if figureIndex==1
    rowIndex = taskId;
else
    rowIndex = taskId-cumulative(figureIndex-1);
end
end

function [inputs,outputs,gains,iterations] = ...
        local_scan(stepMatrices,horizons,stepMs,seed)
inputs = cell(size(horizons));
outputs = cell(size(horizons));
gains = zeros(size(horizons));
iterations = zeros(size(horizons));
currentSeed = seed;
for index = 1:numel(horizons)
    stepCount = round(horizons(index)/stepMs);
    [inputs{index},outputs{index},gains(index),iterations(index)] = ...
        local_top_product_singular(stepMatrices,stepCount,currentSeed);
    if real(inputs{index}'*seed)<0
        inputs{index} = -inputs{index};
        outputs{index} = -outputs{index};
    end
    currentSeed = inputs{index};
end
end

function [right,left,sigma,iteration] = ...
        local_top_product_singular(stepMatrices,stepCount,seed)
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
