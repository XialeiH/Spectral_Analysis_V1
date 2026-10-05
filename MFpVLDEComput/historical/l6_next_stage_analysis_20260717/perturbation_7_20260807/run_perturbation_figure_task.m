function output = run_perturbation_figure_task(caseId, setupFile, outputRoot, ...
        requestedHCOverride,figureNumberOverride,directionSeedOverride)
% Nonlinear ODE perturbation experiments for report Figures 7.1-7.7.

arguments
    caseId (1,1) double {mustBeMember(caseId,1:6)}
    setupFile (1,:) char
    outputRoot (1,:) char
    requestedHCOverride (1,1) double = NaN
    figureNumberOverride (1,:) char = ''
    directionSeedOverride (1,1) double = NaN
end

loaded = load(setupFile, 'setup');
setup = loaded.setup;
context = setup.Context;
fixed = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
tauMs = setup.Config.TauMs;

if caseId == 1
    figureNumber = '7.1';
    l6Weight = 0;
    finalTimeMs = 500;
    sampleStepMs = 2;
    caseLabel = 'Baseline dynamic L6';
    requestedInitialHC = 0.10;
elseif caseId == 2
    figureNumber = '7.2';
    l6Weight = -0.15;
    finalTimeMs = 2500;
    sampleStepMs = 5;
    caseLabel = 'Enhanced L6 derivative gain (1.15x)';
    requestedInitialHC = 0.10;
else
    l6Weight = 0;
    finalTimeMs = 500;
    sampleStepMs = 1;
    if caseId==3
        figureNumber = '7.4';
        caseLabel = 'Baseline dynamic L6, larger structured perturbation';
        requestedInitialHC = 0.30;
        perturbationKind = 'structured';
    else
        figureNumber = sprintf('7.%d',caseId+1);
        requestedInitialHC = [1 5 10];
        requestedInitialHC = requestedInitialHC(caseId-3);
        caseLabel = sprintf('Baseline dynamic L6, random positive perturbation HC %.0f', ...
            requestedInitialHC);
        perturbationKind = 'random_positive';
    end
end
if caseId<=2
    perturbationKind = 'structured';
end
isFixedDirectionOverride = ~isnan(requestedHCOverride);
if isFixedDirectionOverride
    requestedInitialHC = requestedHCOverride;
    figureNumber = figureNumberOverride;
    caseLabel = sprintf(['Baseline dynamic L6, fixed positive direction ' ...
        'perturbation HC %.2f'],requestedInitialHC);
    perturbationKind = 'random_positive';
end

if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

jFixed = l6ns_state_jacobian(setup, 'L6', fixed, l6Weight);
[fixedCluster, fixedEigenvalues, leadingEigenvalue] = local_eigencluster(jFixed);
[fixedSingularOutput, fixedSingularValue, fixedSingularInput] = ...
    svds(jFixed, 1, 'largest');
fixedSingularInput = local_orient_vector(fixedSingularInput, n, context.CWeight);
fixedSingularOutput = local_match_sign(fixedSingularOutput, fixedSingularInput);
maxRealLambda = real(leadingEigenvalue);
if maxRealLambda >= 1
    error('Perturbation7:Unstable', ...
        'Case %d is not stable: max Re lambda = %.9f.', caseId, maxRealLambda);
end

eigenSeedRatio = 0.50;
if strcmp(perturbationKind,'structured')
    eigenSeedDirection = local_orient_vector(real(fixedCluster(:,1)),n,context.CWeight);
    eigenSeedDirection = eigenSeedDirection/max(norm(eigenSeedDirection),eps);
    initialDirection = fixedSingularInput + eigenSeedRatio*eigenSeedDirection;
    initialDirection = initialDirection/max(norm(initialDirection),eps);
else
    directionSeed = 7000+caseId;
    if ~isnan(directionSeedOverride)
        directionSeed = directionSeedOverride;
    end
    rng(directionSeed,'twister');
    initialDirection = abs(randn(size(fixed)));
    initialDirection = initialDirection/max(norm(initialDirection),eps);
    eigenSeedRatio = NaN;
end
unitHC = local_hcnorm(initialDirection, zeros(size(fixed)), n, context.CWeight);
amplitude = requestedInitialHC / max(unitHC, eps);
if strcmp(perturbationKind,'structured')
    negative = initialDirection < 0;
    if any(negative)
        positivityLimit = 0.80 * min(fixed(negative) ./ (-initialDirection(negative)));
        amplitude = min(amplitude, positivityLimit);
    end
end
initialState = fixed + amplitude * initialDirection;
initialHC = local_hcnorm(initialState, fixed, n, context.CWeight);

phi = @(state) l6ns_phi(state, l6Weight, context, [1 1], 1, 'fpp');
rhs = @(~,state) (phi(state)-state) / tauMs;
options = odeset('RelTol',1e-6,'AbsTol',1e-8, ...
    'MaxStep',min(2,sampleStepMs));
fprintf(['Figure %s: w6=%+.3f, max Re lambda=%.9f, sigma1=%.9f, ' ...
    'initial HC=%.6f, kind=%s, T=%.1f ms.\n'], figureNumber, l6Weight, ...
    maxRealLambda, fixedSingularValue, initialHC, perturbationKind, finalTimeMs);
solution = ode45(rhs, [0 finalTimeMs], initialState, options);
times = 0:sampleStepMs:finalTimeMs;
states = deval(solution, times);

delta = states-fixed;
deltaNorm = vecnorm(delta,2,1);
hc = zeros(size(times));
singularAlignment = zeros(size(times));
returnAlignment = zeros(size(times));
for index = 1:numel(times)
    hc(index) = local_hcnorm(states(:,index),fixed,n,context.CWeight);
    denominator = max(deltaNorm(index),eps);
    singularAlignment(index) = abs(fixedSingularOutput' * delta(:,index)) / denominator;
    returnDirection = fixed-states(:,index);
    velocity = phi(states(:,index))-states(:,index);
    returnAlignment(index) = real(returnDirection'*velocity) / ...
        max(norm(returnDirection)*norm(velocity),eps);
end

[transientIndex, convergenceIndex] = local_select_times(times, hc, ...
    deltaNorm, singularAlignment, returnAlignment, tauMs, maxRealLambda);
selectedIndices = [transientIndex convergenceIndex];
selectedTimes = times(selectedIndices);

if caseId>=3
    additionalTimes = [3 75 150];
else
    additionalTimes = [75 150];
end
additionalIndices = zeros(size(additionalTimes));
for index = 1:numel(additionalTimes)
    [~,additionalIndices(index)] = min(abs(times-additionalTimes(index)));
end
allSelectedIndices = unique([selectedIndices additionalIndices],'sorted');
allSnapshots = repmat(struct(),numel(allSelectedIndices),1);
for snapshotIndex = 1:numel(allSelectedIndices)
    stateIndex = allSelectedIndices(snapshotIndex);
    state = states(:,stateIndex);
    jacobian = l6ns_state_jacobian(setup,'L6',state,l6Weight);
    [~,eigenvalues,leading] = local_eigencluster(jacobian);
    [singularOutput,singularValue,singularInput] = svds(jacobian,1,'largest');
    singularInput = local_orient_vector(singularInput,n,context.CWeight);
    singularOutput = local_match_sign(singularOutput,singularInput);
    returnDirection = fixed-state;
    allSnapshots(snapshotIndex).TimeMs = times(stateIndex);
    allSnapshots(snapshotIndex).State = state;
    allSnapshots(snapshotIndex).HCnorm = hc(stateIndex);
    allSnapshots(snapshotIndex).DeltaNorm = deltaNorm(stateIndex);
    allSnapshots(snapshotIndex).SingularAlignment = singularAlignment(stateIndex);
    allSnapshots(snapshotIndex).ReturnAlignment = returnAlignment(stateIndex);
    allSnapshots(snapshotIndex).JacobianMaxReal = real(leading);
    allSnapshots(snapshotIndex).Eigenvalues = eigenvalues;
    allSnapshots(snapshotIndex).ReturnDirectionMaps = ...
        local_signed_maps(returnDirection,mapSize,context.CWeight);
    allSnapshots(snapshotIndex).SingularInputMaps = ...
        local_signed_maps(singularInput,mapSize,context.CWeight);
    allSnapshots(snapshotIndex).SingularOutputMaps = ...
        local_signed_maps(singularOutput,mapSize,context.CWeight);
    allSnapshots(snapshotIndex).FiringMaps = local_state_maps(state,mapSize,context.CWeight);
    allSnapshots(snapshotIndex).TopSingularValue = singularValue;
end
transientSnapshotIndex = find(allSelectedIndices==transientIndex,1);
convergenceSnapshotIndex = find(allSelectedIndices==convergenceIndex,1);
snapshots = allSnapshots([transientSnapshotIndex convergenceSnapshotIndex]);

if isFixedDirectionOverride
    caseTag = sprintf('Baseline_L6_FixedDirection_HC%s', ...
        local_number_tag(requestedInitialHC));
else
    caseTag = local_case_tag(caseId);
end
stem = sprintf('%s_Perturbation_Transient_and_Return_%s',figureNumber,caseTag);

summary = table(string(figureNumber),caseId,string(caseLabel), ...
    string(perturbationKind),l6Weight,1-l6Weight,maxRealLambda, ...
    fixedSingularValue,tauMs,requestedInitialHC,initialHC,eigenSeedRatio, ...
    selectedTimes(1),selectedTimes(2),snapshots(1).HCnorm,snapshots(2).HCnorm, ...
    snapshots(1).SingularAlignment,snapshots(2).SingularAlignment, ...
    snapshots(1).ReturnAlignment,snapshots(2).ReturnAlignment, ...
    snapshots(1).JacobianMaxReal,snapshots(2).JacobianMaxReal, ...
    'VariableNames',{'figureNumber','caseId','caseLabel','perturbationKind', ...
    'l6Weight','l6DerivativeGain','fixedPointMaxRealLambda', ...
    'fixedPointTopSingularValue','tauMs','requestedInitialHC','actualInitialHC', ...
    'eigenSeedRatio','transientTimeMs','convergenceTimeMs','transientHC', ...
    'convergenceHC','transientSingularAlignment','convergenceSingularAlignment', ...
    'transientReturnAlignment','convergenceReturnAlignment', ...
    'transientLocalMaxRealLambda','convergenceLocalMaxRealLambda'});
writetable(summary,fullfile(outputRoot,[stem '_summary.tsv']), ...
    'FileType','text','Delimiter','\t');

snapshotNames = strings(numel(allSnapshots),1);
returnsToFixedPoint = snapshots(2).ReturnAlignment>0 && ...
    snapshots(2).HCnorm<initialHC;
for index = 1:numel(allSnapshots)
    if abs(allSnapshots(index).TimeMs-selectedTimes(1))<1e-9
        snapshotNames(index) = "Ongoing transient";
    elseif abs(allSnapshots(index).TimeMs-selectedTimes(2))<1e-9
        if returnsToFixedPoint
            snapshotNames(index) = "Early return toward fixed point";
        else
            snapshotNames(index) = "Late non-returning state";
        end
    else
        snapshotNames(index) = "Additional snapshot";
    end
end
selectionTable = table((1:numel(allSnapshots))',snapshotNames, ...
    [allSnapshots.TimeMs]',[allSnapshots.HCnorm]', ...
    [allSnapshots.SingularAlignment]',[allSnapshots.ReturnAlignment]', ...
    [allSnapshots.JacobianMaxReal]',[allSnapshots.TopSingularValue]', ...
    'VariableNames',{'snapshot','phase','timeMs','HCnorm', ...
    'fixedSingularOutputAlignment','returnDirectionAlignment', ...
    'localMaxRealLambda','localTopSingularValue'});
writetable(selectionTable,fullfile(outputRoot,[stem '_snapshots.tsv']), ...
    'FileType','text','Delimiter','\t');

output = struct('Summary',summary,'Times',times,'HCnorm',hc, ...
    'DeltaNorm',deltaNorm,'SingularAlignment',singularAlignment, ...
    'ReturnAlignment',returnAlignment,'Snapshots',snapshots, ...
    'AllSnapshots',allSnapshots,'AdditionalSnapshotTimesMs',additionalTimes, ...
    'FixedPointEigenvalues',fixedEigenvalues, ...
    'FixedPointLeadingEigenvalue',leadingEigenvalue, ...
    'FixedPointTopSingularValue',fixedSingularValue, ...
    'InitialState',initialState,'InitialHCnorm',initialHC, ...
    'EigenSeedRatio',eigenSeedRatio,'PerturbationKind',perturbationKind, ...
    'ReturnsToFixedPoint',returnsToFixedPoint, ...
    'L6Weight',l6Weight,'TauMs',tauMs);
save(fullfile(outputRoot,[stem '_data.mat']),'output','-v7.3');
end

function [cluster,eigenvalues,leading] = local_eigencluster(jacobian)
options = struct('tol',1e-9,'maxit',2000,'p',80,'isreal',true);
try
    [vectors,values] = eigs(jacobian,16,'largestreal',options);
catch
    [vectors,values] = eigs(jacobian,16,'lr',options);
end
eigenvalues = diag(values);
[~,order] = sort(real(eigenvalues),'descend');
eigenvalues = eigenvalues(order);
vectors = vectors(:,order);
leading = eigenvalues(1);
clusterTolerance = 5e-4;
clusterCount = sum(real(eigenvalues) >= real(leading)-clusterTolerance);
clusterCount = max(1,min(clusterCount,size(vectors,2)));
cluster = vectors(:,1:clusterCount);
for index = 1:clusterCount
    cluster(:,index) = cluster(:,index)/max(norm(cluster(:,index)),eps);
end
end

function [transientIndex,convergenceIndex] = local_select_times(times,hc, ...
        deltaNorm,singularAlignment,returnAlignment,tauMs,maxRealLambda)
transientCandidates = find(times>=0.5*tauMs & times<=6*tauMs);
if isempty(transientCandidates)
    transientCandidates = 2:min(numel(times),10);
end
[~,localIndex] = max(singularAlignment(transientCandidates));
transientIndex = transientCandidates(localIndex);

slowTimeMs = tauMs/max(1-maxRealLambda,1e-6);
minimumConvergenceTime = times(transientIndex)+max(6*tauMs,0.5*slowTimeMs);
lateCandidates = find(times>=minimumConvergenceTime & ...
    hc>=max(0.02*hc(1),1e-4));
if isempty(lateCandidates)
    lateCandidates = find(times>times(transientIndex));
end
isDecaying = [false diff(deltaNorm)<0];
decayingCandidates = lateCandidates(isDecaying(lateCandidates));
if ~isempty(decayingCandidates)
    lateCandidates = decayingCandidates;
end
lateMaximum = max(returnAlignment(lateCandidates));
qualified = lateCandidates(returnAlignment(lateCandidates)>=0.90*lateMaximum);
if isempty(qualified)
    [~,localIndex] = max(returnAlignment(lateCandidates));
    convergenceIndex = lateCandidates(localIndex);
else
    convergenceIndex = qualified(1);
end
end

function maps = local_state_maps(state,mapSize,wC)
n = prod(mapSize);
s = reshape(state(1:n),mapSize);
c = reshape(state(n+(1:n)),mapSize);
i = reshape(state(2*n+(1:n)),mapSize);
e = (1-wC)*s+wC*c;
maps = {s;c;i;e};
end

function maps = local_signed_maps(vector,mapSize,wC)
maps = local_state_maps(real(vector),mapSize,wC);
for index = 1:numel(maps)
    maps{index} = maps{index}/max(norm(maps{index}(:),2),eps);
end
end


function vector = local_orient_vector(vector,n,wC)
e = (1-wC)*real(vector(1:n))+wC*real(vector(n+(1:n)));
[~,index] = max(abs(e));
if e(index)<0
    vector = -vector;
end
end

function left = local_match_sign(left,right)
if real(left'*right)<0
    left = -left;
end
end

function value = local_hcnorm(state,baseline,n,wC)
dS = state(1:n)-baseline(1:n);
dC = state(n+(1:n))-baseline(n+(1:n));
dI = state(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-wC)*dS+wC*dC;
value = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end

function tag = local_case_tag(caseId)
tags = {'Baseline_L6', 'Enhanced_L6_1p15x', ...
    'Baseline_L6_LargePerturbation_HC0p30', ...
    'Baseline_L6_RandomPositive_HC1', ...
    'Baseline_L6_RandomPositive_HC5', ...
    'Baseline_L6_RandomPositive_HC10'};
tag = tags{caseId};
end

function tag = local_number_tag(value)
tag = strrep(sprintf('%.2f',value),'.','p');
end
