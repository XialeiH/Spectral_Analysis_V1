function result = run_hc10_long_horizon(setupFile, sourceDataFile, outputRoot)
% Continue the deterministic random HC=10 perturbation over a long horizon.

arguments
    setupFile (1,:) char
    sourceDataFile (1,:) char
    outputRoot (1,:) char
end

loadedSetup = load(setupFile,'setup');
setup = loadedSetup.setup;
loadedSource = load(sourceDataFile,'output');
source = loadedSource.output;
context = setup.Context;
fixed = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
tauMs = setup.Config.TauMs;
l6Weight = source.L6Weight;
phi = @(state) l6ns_phi(state,l6Weight,context,[1 1],1,'fpp');
rhs = @(~,state) (phi(state)-state)/tauMs;

checkpointTimesMs = [0 100 500 1000 2000 5000 10000 20000];
states = zeros(numel(fixed),numel(checkpointTimesMs));
states(:,1) = source.InitialState;
options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',10);
for index = 2:numel(checkpointTimesMs)
    interval = checkpointTimesMs(index-1:index);
    fprintf('Integrating HC=10 trajectory from %.0f to %.0f ms.\n',interval);
    solution = ode45(rhs,interval,states(:,index-1),options);
    states(:,index) = deval(solution,interval(2));
end

hcFromOriginal = zeros(size(checkpointTimesMs));
fixedPointResidualHC = zeros(size(checkpointTimesMs));
returnAlignment = zeros(size(checkpointTimesMs));
stateChangeHC = NaN(size(checkpointTimesMs));
minimumRate = zeros(size(checkpointTimesMs));
maximumRate = zeros(size(checkpointTimesMs));
for index = 1:numel(checkpointTimesMs)
    state = states(:,index);
    mapped = phi(state);
    velocityNumerator = mapped-state;
    returnDirection = fixed-state;
    hcFromOriginal(index) = local_hcnorm(state,fixed,n,context.CWeight);
    fixedPointResidualHC(index) = local_hcnorm(mapped,state,n,context.CWeight);
    returnAlignment(index) = real(returnDirection'*velocityNumerator)/ ...
        max(norm(returnDirection)*norm(velocityNumerator),eps);
    if index>1
        stateChangeHC(index) = local_hcnorm(state,states(:,index-1),n,context.CWeight);
    end
    minimumRate(index) = min(state);
    maximumRate(index) = max(state);
end

iterationState = states(:,end);
iterationResidualHC = zeros(100,1);
for index = 1:100
    mapped = phi(iterationState);
    iterationResidualHC(index) = local_hcnorm(mapped,iterationState,n,context.CWeight);
    iterationState = mapped;
end

summary = table(checkpointTimesMs(:),hcFromOriginal(:), ...
    fixedPointResidualHC(:),returnAlignment(:),stateChangeHC(:), ...
    minimumRate(:),maximumRate(:), ...
    'VariableNames',{'timeMs','HCFromOriginalFixedPoint', ...
    'fixedPointResidualHC','returnAlignment','HCChangeFromPreviousCheckpoint', ...
    'minimumRate','maximumRate'});

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
writetable(summary,fullfile(outputRoot,'hc10_long_horizon_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
result = struct('Summary',summary,'States',states, ...
    'IterationResidualHC',iterationResidualHC, ...
    'IterationFinalHCFromOriginal', ...
    local_hcnorm(iterationState,fixed,n,context.CWeight), ...
    'IterationFinalHCFromODE20s', ...
    local_hcnorm(iterationState,states(:,end),n,context.CWeight), ...
    'TauMs',tauMs,'L6Weight',l6Weight,'MapSize',mapSize);
save(fullfile(outputRoot,'hc10_long_horizon_data.mat'),'result','-v7.3');
fprintf(['At 20 s: HC(original)=%.9g, residual HC=%.9g; ' ...
    'after 100 map iterations: HC(original)=%.9g, HC(from ODE20s)=%.9g.\n'], ...
    hcFromOriginal(end),fixedPointResidualHC(end), ...
    result.IterationFinalHCFromOriginal,result.IterationFinalHCFromODE20s);
end

function value = local_hcnorm(state,baseline,n,wC)
dS = state(1:n)-baseline(1:n);
dC = state(n+(1:n))-baseline(n+(1:n));
dI = state(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-wC)*dS+wC*dC;
value = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end
