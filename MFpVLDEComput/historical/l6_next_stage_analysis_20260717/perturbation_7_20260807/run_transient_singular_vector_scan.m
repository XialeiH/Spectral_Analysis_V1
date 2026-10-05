function result = run_transient_singular_vector_scan(caseId,setupFile,outputRoot)
% Track the local top right singular vector along the early ODE transient.

arguments
    caseId (1,1) double {mustBeMember(caseId,[1 2])}
    setupFile (1,:) char
    outputRoot (1,:) char
end

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
fixed = context.FixedPoint(:);
mapSize = double(context.MapSize(:).');
n = prod(mapSize);
tauMs = setup.Config.TauMs;

if caseId==1
    l6Weight = 0;
    caseLabel = 'Baseline dynamic L6';
else
    l6Weight = -0.15;
    caseLabel = 'Enhanced L6 derivative gain (1.15x)';
end

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end

jFixed = l6ns_state_jacobian(setup,'L6',fixed,l6Weight);
[fixedCluster,~,~] = local_eigencluster(jFixed);
[~,~,fixedSingularInput] = svds(jFixed,1,'largest');
fixedSingularInput = local_orient_vector(fixedSingularInput,n,context.CWeight);

requestedInitialHC = 0.10;
eigenSeedRatio = 0.50;
eigenSeedDirection = local_orient_vector(real(fixedCluster(:,1)),n,context.CWeight);
eigenSeedDirection = eigenSeedDirection/max(norm(eigenSeedDirection),eps);
initialDirection = fixedSingularInput + eigenSeedRatio*eigenSeedDirection;
initialDirection = initialDirection/max(norm(initialDirection),eps);
unitHC = local_hcnorm(initialDirection,zeros(size(fixed)),n,context.CWeight);
amplitude = requestedInitialHC/max(unitHC,eps);
negative = initialDirection<0;
if any(negative)
    positivityLimit = 0.80*min(fixed(negative)./(-initialDirection(negative)));
    amplitude = min(amplitude,positivityLimit);
end
initialState = fixed + amplitude*initialDirection;

scanTimesMs = 0:5:250;
phi = @(state) l6ns_phi(state,l6Weight,context,[1 1],1,'fpp');
rhs = @(~,state) (phi(state)-state)/tauMs;
options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',2);
solution = ode45(rhs,[0 scanTimesMs(end)],initialState,options);
states = deval(solution,scanTimesMs);

count = numel(scanTimesMs);
sigma1 = zeros(count,1);
sigma2 = zeros(count,1);
gapRatio = zeros(count,1);
referenceCosine = zeros(count,1);
consecutiveCosine = ones(count,1);
referenceEMapCosine = zeros(count,1);
consecutiveEMapCosine = ones(count,1);
referenceRelativeChange = zeros(count,1);
consecutiveRelativeChange = zeros(count,1);
top2SubspaceOverlap = zeros(count,1);
eMaps = zeros([mapSize count]);

referenceVector = [];
referenceBasis = [];
referenceE = [];
previousVector = [];
previousE = [];

for index = 1:count
    fprintf('Case %d singular scan %d/%d, t=%.1f ms.\n', ...
        caseId,index,count,scanTimesMs(index));
    jacobian = l6ns_state_jacobian(setup,'L6',states(:,index),l6Weight);
    [~,singularValues,rightVectors] = svds(jacobian,2,'largest');
    singularValues = diag(singularValues);
    [singularValues,order] = sort(real(singularValues),'descend');
    rightVectors = real(rightVectors(:,order));
    for modeIndex = 1:2
        rightVectors(:,modeIndex) = rightVectors(:,modeIndex) / ...
            max(norm(rightVectors(:,modeIndex)),eps);
    end
    vector = local_orient_vector(rightVectors(:,1),n,context.CWeight);
    if dot(vector,rightVectors(:,1))<0
        rightVectors(:,1) = -rightVectors(:,1);
    end
    eVector = local_e_vector(vector,n,context.CWeight);
    eVector = eVector/max(norm(eVector),eps);

    if index==1
        referenceVector = vector;
        referenceBasis = orth(rightVectors(:,1:2));
        referenceE = eVector;
    end
    if index>1
        if dot(vector,previousVector)<0
            vector = -vector;
            eVector = -eVector;
        end
        consecutiveCosine(index) = dot(previousVector,vector);
        consecutiveEMapCosine(index) = dot(previousE,eVector);
        consecutiveRelativeChange(index) = norm(eVector-previousE) / ...
            max(norm(previousE),eps);
    end

    if dot(vector,referenceVector)<0
        referenceAlignedVector = -vector;
        referenceAlignedE = -eVector;
    else
        referenceAlignedVector = vector;
        referenceAlignedE = eVector;
    end
    currentBasis = orth(rightVectors(:,1:2));
    sigma1(index) = singularValues(1);
    sigma2(index) = singularValues(2);
    gapRatio(index) = (sigma1(index)-sigma2(index))/max(sigma1(index),eps);
    referenceCosine(index) = dot(referenceVector,referenceAlignedVector);
    referenceEMapCosine(index) = dot(referenceE,referenceAlignedE);
    referenceRelativeChange(index) = norm(referenceAlignedE-referenceE) / ...
        max(norm(referenceE),eps);
    top2SubspaceOverlap(index) = norm(referenceBasis'*currentBasis,'fro')/sqrt(2);
    eMaps(:,:,index) = reshape(referenceAlignedE,mapSize);
    previousVector = vector;
    previousE = eVector;
end

metrics = table(scanTimesMs(:),sigma1,sigma2,gapRatio,referenceCosine, ...
    consecutiveCosine,referenceEMapCosine,consecutiveEMapCosine, ...
    referenceRelativeChange,consecutiveRelativeChange,top2SubspaceOverlap, ...
    'VariableNames',{'timeMs','sigma1','sigma2','gapRatio', ...
    'fullVectorCosineToT0','fullVectorConsecutiveCosine', ...
    'eMapCosineToT0','eMapConsecutiveCosine','eMapRelativeChangeToT0', ...
    'eMapConsecutiveRelativeChange','top2SubspaceOverlapToT0'});

result = struct('CaseId',caseId,'CaseLabel',caseLabel,'L6Weight',l6Weight, ...
    'TauMs',tauMs,'MapSize',mapSize,'TimesMs',scanTimesMs, ...
    'Metrics',metrics,'ESingularMaps',eMaps);
stem = sprintf('transient_singular_scan_case%d',caseId);
save(fullfile(outputRoot,[stem '.mat']),'result','-v7.3');
writetable(metrics,fullfile(outputRoot,[stem '.tsv']), ...
    'FileType','text','Delimiter','\t');
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
clusterCount = sum(real(eigenvalues)>=real(leading)-5e-4);
clusterCount = max(1,min(clusterCount,size(vectors,2)));
cluster = vectors(:,1:clusterCount);
end

function vector = local_orient_vector(vector,n,wC)
vector = real(vector)/max(norm(vector),eps);
e = local_e_vector(vector,n,wC);
[~,index] = max(abs(e));
if e(index)<0
    vector = -vector;
end
end

function e = local_e_vector(vector,n,wC)
e = (1-wC)*real(vector(1:n)) + wC*real(vector(n+(1:n)));
end

function value = local_hcnorm(state,baseline,n,wC)
dS = state(1:n)-baseline(1:n);
dC = state(n+(1:n))-baseline(n+(1:n));
dI = state(2*n+(1:n))-baseline(2*n+(1:n));
dE = (1-wC)*dS+wC*dC;
value = sqrt(mean(0.8*dE.^2+0.2*dI.^2));
end
