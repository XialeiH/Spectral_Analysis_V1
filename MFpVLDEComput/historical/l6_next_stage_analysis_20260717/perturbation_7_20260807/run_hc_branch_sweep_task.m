function result = run_hc_branch_sweep_task(taskId,setupFile,outputRoot)
% Sweep perturbation HC norm along one fixed positive direction.

arguments
    taskId (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    outputRoot (1,:) char
end

hcValues = [0.10 0.30 0.50 0.75 1.00 1.50 2.00 2.50 3.00 3.50 ...
    4.00 4.50 5.00 5.50 6.00 6.50 7.00 7.50 8.00 9.00 10.00 ...
    6.10 6.20 6.25 6.30 6.35 6.40 6.45 ...
    6.01 6.02 6.03 6.04 6.05 6.06 6.07 6.08 6.09];
if taskId>numel(hcValues)
    error('Perturbation7:TaskId','Task %d exceeds the %d HC values.', ...
        taskId,numel(hcValues));
end
requestedHC = hcValues(taskId);

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
fixed = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
tauMs = setup.Config.TauMs;

% Use the former HC=10 experiment's direction for every amplitude.
rng(7006,'twister');
direction = abs(randn(size(fixed)));
direction = direction/max(norm(direction),eps);
unitHC = local_hcnorm(direction,zeros(size(fixed)),n,context.CWeight);
initialState = fixed+(requestedHC/unitHC)*direction;
initialHC = local_hcnorm(initialState,fixed,n,context.CWeight);

phi = @(state) l6ns_phi(state,0,context,[1 1],1,'fpp');
rhs = @(~,state) (phi(state)-state)/tauMs;
times = unique([0:2:500 510:10:1000 1025:25:2000]);
options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',10);
fprintf('Fixed-direction HC sweep task %d: HC %.3f to %.0f ms.\n', ...
    taskId,initialHC,times(end));
solution = ode45(rhs,[times(1) times(end)],initialState,options);
states = deval(solution,times);

hcTrajectory = zeros(size(times));
minimumRate = zeros(size(times));
maximumRate = zeros(size(times));
for index = 1:numel(times)
    hcTrajectory(index) = local_hcnorm(states(:,index),fixed,n,context.CWeight);
    minimumRate(index) = min(states(:,index));
    maximumRate(index) = max(states(:,index));
end
terminalState = states(:,end);
terminalMapped = phi(terminalState);
terminalResidualHC = local_hcnorm(terminalMapped,terminalState,n,context.CWeight);
terminalHC = hcTrajectory(end);
terminalL4E = local_l4e_inputs(terminalState,context);
l4eMinimum = min(cellfun(@(value) min(value(:)),terminalL4E));
l4eMaximum = max(cellfun(@(value) max(value(:)),terminalL4E));
fractionAbove47500 = mean(vertcat(terminalL4E{:})>47500);
fractionAbove55000 = mean(vertcat(terminalL4E{:})>55000);

if terminalResidualHC<1e-3 && fractionAbove55000>0.99
    regime = "artificial high fixed point";
elseif terminalResidualHC<1e-3 && terminalHC<1
    regime = "physiological fixed point";
else
    regime = "not converged by 2000 ms";
end

initialE = local_e_map(initialState,mapSize,context.CWeight);
terminalE = local_e_map(terminalState,mapSize,context.CWeight);
result = struct('TaskId',taskId,'RequestedHC',requestedHC, ...
    'InitialHC',initialHC,'TimesMs',times,'HCtrajectory',hcTrajectory, ...
    'MinimumRate',minimumRate,'MaximumRate',maximumRate, ...
    'TerminalHCFromOriginal',terminalHC, ...
    'TerminalFixedPointResidualHC',terminalResidualHC, ...
    'TerminalL4EMinimum',l4eMinimum,'TerminalL4EMaximum',l4eMaximum, ...
    'FractionL4EAbove47500',fractionAbove47500, ...
    'FractionL4EAbove55000',fractionAbove55000,'Regime',regime, ...
    'InitialEMap',initialE,'TerminalEMap',terminalE,'TauMs',tauMs, ...
    'DirectionSeed',7006,'FinalTimeMs',times(end));
if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
save(fullfile(outputRoot,sprintf('hc_branch_sweep_%02d.mat',taskId)), ...
    'result','-v7.3');
fprintf(['HC %.3f: terminal HC %.6g, residual %.6g, L4E min %.6g, ' ...
    'above55000 %.3f, %s.\n'],initialHC,terminalHC,terminalResidualHC, ...
    l4eMinimum,fractionAbove55000,regime);
end

function values = local_l4e_inputs(state,context)
n = numel(state)/3;
s = state(1:n);
c = state(n+(1:n));
i = state(2*n+(1:n));
if context.Isaturation
    eRaw = (1-context.CWeight)*s+context.CWeight*c;
    eBase = L6Convert(eRaw,context.EKpUse);
    adjustment = eBase./eRaw;
    adjustment(~isfinite(adjustment)) = 1;
    s = s.*adjustment;
    c = c.*adjustment;
    i = InhMulp(i,context.IKpUse);
end
values = {(context.C_SS*s+context.C_SC*c)/context.L4SEp; ...
    (context.C_CS*s+context.C_CC*c)/context.L4CEp; ...
    (context.C_IS*s+context.C_IC*c)/context.L4IEp};
end

function map = local_e_map(state,mapSize,wC)
n = prod(mapSize);
s = reshape(state(1:n),mapSize);
c = reshape(state(n+(1:n)),mapSize);
map = (1-wC)*s+wC*c;
end

function value = local_hcnorm(state,baseline,n,wC)
dS = state(1:n)-baseline(1:n);
dC = state(n+(1:n))-baseline(n+(1:n));
dI = state(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-wC)*dS+wC*dC;
value = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end
