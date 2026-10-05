function result = compare_hc5_snapshot_finite_time_transient(setupFile,outputRoot)
% Compare snapshot-Jacobian and trajectory-aware finite-time transient modes.

arguments
    setupFile (1,:) char
    outputRoot (1,:) char
end

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
fixed = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
dimension = numel(fixed);
tauMs = setup.Config.TauMs;
l6Weight = 0;

% Reproduce the existing HC=5 random-positive perturbation exactly.
rng(7005,'twister');
initialDirection = abs(randn(size(fixed)));
initialDirection = initialDirection/max(norm(initialDirection),eps);
unitHC = local_hcnorm(initialDirection,zeros(size(fixed)),n,context.CWeight);
initialState = fixed + (5/unitHC)*initialDirection;
initialHC = local_hcnorm(initialState,fixed,n,context.CWeight);

phi = @(state) l6ns_phi(state,l6Weight,context,[1 1],1,'fpp');
rhs = @(~,state) (phi(state)-state)/tauMs;
horizonGrid = [0.5 1 2 3 5 8 12 15 20 30 45 60 80 100];
stepMs = 0.5;
timeGrid = 0:stepMs:max(horizonGrid);
options = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);
solution = ode45(rhs,[0 max(timeGrid)],initialState,options);
states = deval(solution,timeGrid);

% Current method: SVD of the instantaneous D Phi at the HC=5 state.
jSnapshot = sparse(l6ns_state_jacobian(setup,'L6',initialState,l6Weight));
[snapshotOutput,snapshotSigma,snapshotInput] = svds(jSnapshot,1,'largest');
snapshotInput = local_orient_e(snapshotInput,n,context.CWeight);
if real(snapshotOutput'*snapshotInput)<0
    snapshotOutput = -snapshotOutput;
end

% Proposed method: time-ordered tangent propagator for
% A(t)=(D Phi(x(t))-I)/tau. Forward Euler gives an exactly adjointable
% discrete product and stepMs/tau is below 0.05 here.
identity = speye(dimension);
stepMatrices = cell(numel(timeGrid)-1,1);
for index = 1:numel(stepMatrices)
    jacobian = sparse(l6ns_state_jacobian( ...
        setup,'L6',states(:,index),l6Weight));
    stepMatrices{index} = identity + ...
        (stepMs/tauMs)*(jacobian-identity);
    if mod(index,20)==0 || index==numel(stepMatrices)
        fprintf('Built tangent step %d/%d at %.1f ms.\n', ...
            index,numel(stepMatrices),timeGrid(index));
    end
end

finiteGain = zeros(size(horizonGrid));
finiteInput = cell(size(horizonGrid));
finiteOutput = cell(size(horizonGrid));
iterations = zeros(size(horizonGrid));
seed = snapshotInput;
for index = 1:numel(horizonGrid)
    stepCount = round(horizonGrid(index)/stepMs);
    [finiteInput{index},finiteOutput{index},finiteGain(index),iterations(index)] = ...
        local_top_product_singular(stepMatrices,stepCount,seed);
    if real(finiteInput{index}'*snapshotInput)<0
        finiteInput{index} = -finiteInput{index};
        finiteOutput{index} = -finiteOutput{index};
    end
    seed = finiteInput{index};
    fprintf('T %.1f ms: finite-time gain %.9g (%d iterations).\n', ...
        horizonGrid(index),finiteGain(index),iterations(index));
end

[peakGain,peakIndex] = max(finiteGain);
peakTimeMs = horizonGrid(peakIndex);
optimalInput = finiteInput{peakIndex};
optimalOutput = finiteOutput{peakIndex};
inputCosine = abs(real(snapshotInput'*optimalInput)) / ...
    max(norm(snapshotInput)*norm(optimalInput),eps);
outputCosine = abs(real(snapshotOutput'*optimalOutput)) / ...
    max(norm(snapshotOutput)*norm(optimalOutput),eps);

snapshotInputE = local_e_map(snapshotInput,mapSize,context.CWeight);
optimalInputE = local_e_map(optimalInput,mapSize,context.CWeight);
snapshotOutputE = local_e_map(snapshotOutput,mapSize,context.CWeight);
optimalOutputE = local_e_map(optimalOutput,mapSize,context.CWeight);
inputDifferenceE = optimalInputE-snapshotInputE;
outputDifferenceE = optimalOutputE-snapshotOutputE;

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
figureHandle = figure('Visible','off','Color','w', ...
    'Position',[80 80 1700 900]);
layout = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

inputLimit = max(abs([snapshotInputE(:);optimalInputE(:)]));
outputLimit = max(abs([snapshotOutputE(:);optimalOutputE(:)]));
inputDiffLimit = max(abs(inputDifferenceE(:)));
outputDiffLimit = max(abs(outputDifferenceE(:)));

local_map(nexttile(layout,1),snapshotInputE,inputLimit, ...
    'Current: top right singular vector of D\Phi(x_0)');
local_map(nexttile(layout,2),optimalInputE,inputLimit, ...
    sprintf('Finite-time optimal input, T_{peak}=%.1f ms',peakTimeMs));
local_map(nexttile(layout,3),inputDifferenceE,inputDiffLimit, ...
    sprintf('Input difference; |cos|=%.3f',inputCosine));
axisGain = nexttile(layout,4);
plot(axisGain,horizonGrid,finiteGain,'-o','LineWidth',1.8, ...
    'MarkerFaceColor',[0.00 0.45 0.74]); hold(axisGain,'on');
xline(axisGain,peakTimeMs,'--r',sprintf('T_{peak}=%.1f ms',peakTimeMs), ...
    'LabelVerticalAlignment','bottom');
xlabel(axisGain,'Finite-time horizon T (ms)');
ylabel(axisGain,'\sigma_{max}(P(T,0))');
title(axisGain,'Trajectory-aware finite-time amplification');
grid(axisGain,'on'); box(axisGain,'on');

local_map(nexttile(layout,5),snapshotOutputE,outputLimit, ...
    sprintf('Current: top left singular vector of D\Phi(x_0), \sigma=%.2f',snapshotSigma));
local_map(nexttile(layout,6),optimalOutputE,outputLimit, ...
    sprintf('Finite-time response at %.1f ms, gain=%.2f',peakTimeMs,peakGain));
local_map(nexttile(layout,7),outputDifferenceE,outputDiffLimit, ...
    sprintf('Output difference; |cos|=%.3f',outputCosine));
axisSummary = nexttile(layout,8);
axis(axisSummary,'off');
text(axisSummary,0.02,0.92,{ ...
    'Comparison definition', ...
    sprintf('Base trajectory: random-positive HC norm = %.1f',initialHC), ...
    'Current: SVD of instantaneous D\Phi(x_0)', ...
    'Proposed: SVD of time-ordered P(T,0)', ...
    sprintf('Peak horizon: %.1f ms',peakTimeMs), ...
    sprintf('Peak finite-time gain: %.4g',peakGain), ...
    sprintf('Input full-state |cosine|: %.4f',inputCosine), ...
    sprintf('Output full-state |cosine|: %.4f',outputCosine), ...
    sprintf('Tangent discretization: %.1f ms',stepMs)}, ...
    'Units','normalized','VerticalAlignment','top','FontSize',11);

title(layout,{ ...
    'HC=5 transient direction: instantaneous Jacobian SVD versus finite-time propagator SVD', ...
    'E population maps; contrast 100, orientation 0.00 deg, baseline dynamic L6'}, ...
    'FontWeight','bold','FontSize',14);
pdfFile = fullfile(outputRoot, ...
    'HC5_Instantaneous_vs_FiniteTime_Transient_Direction.pdf');
exportgraphics(figureHandle,pdfFile,'ContentType','vector');
close(figureHandle);

summary = table(initialHC,snapshotSigma,peakTimeMs,peakGain,inputCosine, ...
    outputCosine,stepMs,string(pdfFile), ...
    'VariableNames',{'initialHCnorm','snapshotDPhiSigma1','peakTimeMs', ...
    'finiteTimePeakGain','inputAbsoluteCosine','outputAbsoluteCosine', ...
    'tangentStepMs','pdfFile'});
writetable(summary,fullfile(outputRoot, ...
    'HC5_Instantaneous_vs_FiniteTime_Transient_Direction_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

result = struct('Summary',summary,'HorizonGridMs',horizonGrid, ...
    'FiniteTimeGain',finiteGain,'PowerIterations',iterations, ...
    'SnapshotInput',snapshotInput,'SnapshotOutput',snapshotOutput, ...
    'FiniteTimeOptimalInput',optimalInput, ...
    'FiniteTimeOptimalOutput',optimalOutput,'PeakTimeMs',peakTimeMs, ...
    'InitialState',initialState,'InitialHCnorm',initialHC, ...
    'TimeGridMs',timeGrid,'TrajectoryHCnorm', ...
    local_trajectory_hc(states,fixed,n,context.CWeight));
save(fullfile(outputRoot, ...
    'HC5_Instantaneous_vs_FiniteTime_Transient_Direction_data.mat'), ...
    'result','-v7.3');
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

function map = local_e_map(vector,mapSize,wC)
n = prod(mapSize);
map = reshape((1-wC)*real(vector(1:n))+ ...
    wC*real(vector(n+(1:n))),mapSize);
end

function vector = local_orient_e(vector,n,wC)
e = (1-wC)*real(vector(1:n))+wC*real(vector(n+(1:n)));
[~,index] = max(abs(e));
if e(index)<0
    vector = -vector;
end
end

function local_map(axisHandle,map,limitValue,titleText)
imagesc(axisHandle,map);
axis(axisHandle,'image');
set(axisHandle,'YDir','normal','XTick',[],'YTick',[]);
colormap(axisHandle,jet(256));
limitValue = max(limitValue,eps);
clim(axisHandle,[-limitValue limitValue]);
colorbar(axisHandle);
title(axisHandle,titleText,'FontSize',10);
end

function values = local_trajectory_hc(states,fixed,n,wC)
values = zeros(1,size(states,2));
for index = 1:size(states,2)
    values(index) = local_hcnorm(states(:,index),fixed,n,wC);
end
end

function value = local_hcnorm(state,baseline,n,wC)
dS = state(1:n)-baseline(1:n);
dC = state(n+(1:n))-baseline(n+(1:n));
dI = state(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-wC)*dS+wC*dC;
value = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end
