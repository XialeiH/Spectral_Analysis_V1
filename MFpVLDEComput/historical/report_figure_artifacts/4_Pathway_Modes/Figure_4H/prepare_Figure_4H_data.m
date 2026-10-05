function prepare_Figure_4H_data()
% Compute decoder trajectories and validation metrics for Figure 4H.

artifactRoot = fileparts(mfilename('fullpath'));
projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project';
projectRootOverride = getenv('FIGURE4H_PROJECT_ROOT');
if ~isempty(projectRootOverride); projectRoot = projectRootOverride; end
modelRoot = fullfile(projectRoot,'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model','Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main');
runtimeRoot = fullfile(projectRoot,'Spectral_Analysis', ...
    'sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134', ...
    'hpc_runtime_bundle');
l6Root = fullfile(projectRoot,'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model','l6_next_stage_analysis_20260717');
modelRootOverride = getenv('FIGURE4H_MODEL_ROOT');
runtimeRootOverride = getenv('FIGURE4H_RUNTIME_ROOT');
l6RootOverride = getenv('FIGURE4H_L6_ROOT');
if ~isempty(modelRootOverride); modelRoot = modelRootOverride; end
if ~isempty(runtimeRootOverride); runtimeRoot = runtimeRootOverride; end
if ~isempty(l6RootOverride); l6Root = l6RootOverride; end
addpath(fullfile(modelRoot,'Utils'),'-begin');
addpath(runtimeRoot,'-begin');
addpath(l6Root,'-begin');

dataRoot = fullfile(artifactRoot,'data');
setupData = load(fullfile(dataRoot,'global_bifurcation_setup.mat'),'setup');
geometry = load(fullfile(dataRoot,'Figure_4H_geometry.mat'), ...
    'qOrientation','qLogContrast');
trajectory = load(fullfile(dataRoot,'figure_04_trajectory.mat'), ...
    'states','trajectoryTimesMs');
setup = setupData.setup;
context = setup.Context;
fixedPoint = double(context.FixedPoint(:));
tauMs = double(setup.Config.TauMs);
mapSize = double(context.MapSize(:).');

% Remove gain and contrast directions before decoding orientation.
qOrientation = double(geometry.qOrientation(:));
qLogContrast = double(geometry.qLogContrast(:));
nuisance = orth([fixedPoint,qLogContrast]);
qTheta = qOrientation-nuisance*(nuisance'*qOrientation);
decoderDenominator = qTheta'*qTheta;
assert(decoderDenominator>eps,'Degenerate orientation decoder.');

timeMask = trajectory.trajectoryTimesMs<=180;
timesMs = double(trajectory.trajectoryTimesMs(timeMask));
dnnStates = double(trajectory.states(:,timeMask));
initialDisplacement = dnnStates(:,1)-fixedPoint;

% Fixed-J linear prediction for the identical initial perturbation.
generator = sparse(double(setup.Pathway.JBaseline)-speye(numel(fixedPoint)))/tauMs;
timeStepMs = timesMs(2)-timesMs(1);
linearStep = speye(numel(fixedPoint))+timeStepMs*generator;
linearDelta = zeros(numel(fixedPoint),numel(timesMs));
linearDelta(:,1) = initialDisplacement;
for timeIndex = 2:numel(timesMs)
    linearDelta(:,timeIndex) = linearStep*linearDelta(:,timeIndex-1);
end

regionNames = ["Optimal","Neighboring","Orthogonal"];
spatialRegion = local_orientation_regions(mapSize);
stateRegions = false(numel(fixedPoint),3);
for region = 1:3
    spatialMask = spatialRegion==region;
    stateRegions(:,region) = repmat(spatialMask(:),3,1);
end

dnnDelta = dnnStates-fixedPoint;
dnnDecoder = local_decode(qTheta,decoderDenominator,dnnDelta);
linearDecoder = local_decode(qTheta,decoderDenominator,linearDelta);
dnnRegional = local_regional_decode(qTheta,decoderDenominator,dnnDelta,stateRegions);
linearRegional = local_regional_decode(qTheta,decoderDenominator,linearDelta,stateRegions);
assert(max(abs(sum(dnnRegional,2)-dnnDecoder))<1e-9, ...
    'Regional DNN contributions do not sum to the full decoder.');

% Parent CG uses the same readout and matched HC=1 displacement.
parentCGFile = fullfile(dataRoot,'Figure_4H_parent_CG.mat');
parentCGAvailable = exist(parentCGFile,'file')==2;
parentCGTimesMs = [];
parentCGDecoder = [];
parentCGRegional = [];
parentCGIterationMs = NaN;
if parentCGAvailable
    parent = load(parentCGFile,'parentCGStates','parentCGTimesMs', ...
        'parentCGFixedVector');
    % The parent library uses x_{n+1}=x_n+p(Phi(x_n)-x_n), p=0.33.
    % With tau calibrated in the continuous ODE, one library epoch is p*tau ms.
    parentCGIterationMs = 0.33*tauMs;
    parentCGTimesMs = double(parent.parentCGTimesMs(:).')*parentCGIterationMs;
    parentMask = parentCGTimesMs<=180;
    parentCGTimesMs = parentCGTimesMs(parentMask);
    parentDelta = double(parent.parentCGStates(:,parentMask))- ...
        double(parent.parentCGFixedVector(:));
    parentCGDecoder = local_decode(qTheta,decoderDenominator,parentDelta);
    parentCGRegional = local_regional_decode( ...
        qTheta,decoderDenominator,parentDelta,stateRegions);
end

% Nonlinear amplitude and sign check on a coarser time grid.
validationTimesMs = 0:2:180;
amplitudeScales = [0.25 0.5 1.0];
validationDecoder = zeros(numel(validationTimesMs),numel(amplitudeScales),2);
phi = @(state)l6ns_phi(state,0,context,[1 1],1,'fpp');
nonlinearOptions = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1);
for scaleIndex = 1:numel(amplitudeScales)
    for signIndex = 1:2
        signValue = [-1 1];
        initialState = fixedPoint+signValue(signIndex)* ...
            amplitudeScales(scaleIndex)*initialDisplacement;
        solution = ode45(@(~,state)(phi(state)-state)/tauMs, ...
            [0 validationTimesMs(end)],initialState,nonlinearOptions);
        states = deval(solution,validationTimesMs);
        validationDecoder(:,scaleIndex,signIndex) = local_decode( ...
            qTheta,decoderDenominator,states-fixedPoint);
    end
end
normalizedPositive = squeeze(validationDecoder(:,:,2))./amplitudeScales;
normalizedNegative = -squeeze(validationDecoder(:,:,1))./amplitudeScales;
linearityReference = normalizedPositive(:,1);
linearityRelativeRms = zeros(size(amplitudeScales));
signSymmetryRelativeRms = zeros(size(amplitudeScales));
for index = 1:numel(amplitudeScales)
    scaleDenominator = max(rms(linearityReference),eps);
    linearityRelativeRms(index) = rms( ...
        normalizedPositive(:,index)-linearityReference)/scaleDenominator;
    signSymmetryRelativeRms(index) = rms( ...
        normalizedPositive(:,index)-normalizedNegative(:,index))/ ...
        max(rms(normalizedPositive(:,index)),eps);
end

dnnMetrics = local_metrics(timesMs,dnnDecoder,dnnRegional);
linearMetrics = local_metrics(timesMs,linearDecoder,linearRegional);
if parentCGAvailable
    parentCGMetrics = local_metrics(parentCGTimesMs,parentCGDecoder,parentCGRegional);
else
    parentCGMetrics = struct();
end

save(fullfile(dataRoot,'Figure_4H_analysis.mat'), ...
    'timesMs','dnnDecoder','linearDecoder','dnnRegional','linearRegional', ...
    'parentCGAvailable','parentCGTimesMs','parentCGDecoder','parentCGRegional', ...
    'parentCGIterationMs', ...
    'regionNames','spatialRegion','qTheta','decoderDenominator','tauMs', ...
    'validationTimesMs','amplitudeScales','validationDecoder', ...
    'linearityRelativeRms','signSymmetryRelativeRms', ...
    'dnnMetrics','linearMetrics','parentCGMetrics','-v7.3');
fprintf('Saved Figure 4H analysis. DNN peak %.4f deg at %.1f ms.\n', ...
    dnnMetrics.PeakErrorDeg,dnnMetrics.PeakTimeMs);
fprintf('Linearity RMS: %s; sign-symmetry RMS: %s.\n', ...
    mat2str(linearityRelativeRms,4),mat2str(signSymmetryRelativeRms,4));
end

function decoded = local_decode(qTheta,denominator,delta)
decoded = (qTheta'*delta/denominator).';
end

function regional = local_regional_decode(qTheta,denominator,delta,masks)
regional = zeros(size(delta,2),size(masks,2));
for region = 1:size(masks,2)
    maskedQ = qTheta.*masks(:,region);
    regional(:,region) = (maskedQ'*delta/denominator).';
end
end

function regions = local_orientation_regions(mapSize)
% Canonical pinwheel geometry, defined before inspecting the perturbation.
[columns,rows] = meshgrid(1:mapSize(2),1:mapSize(1));
localColumn = mod(columns-1,10)+1;
localRow = mod(rows-1,10)+1;
angle = mod(0.5*atan2d(5.5-localRow,localColumn-5.5),180);
distance = min(angle,180-angle);
regions = ones(mapSize);
regions(distance>15 & distance<67.5) = 2;
regions(distance>=67.5) = 3;
end

function metrics = local_metrics(times,decoder,regional)
[peakValue,peakIndex] = max(abs(decoder));
metrics.PeakErrorDeg = peakValue;
metrics.PeakTimeMs = times(peakIndex);
metrics.LateResidualDeg = mean(abs(decoder(times>=120)));
metrics.RegionPeakDeg = zeros(1,3);
metrics.RegionPeakTimeMs = zeros(1,3);
metrics.RegionPersistenceMs = nan(1,3);
metrics.RegionLateResidualDeg = zeros(1,3);
for region = 1:3
    [regionPeak,index] = max(abs(regional(:,region)));
    metrics.RegionPeakDeg(region) = regionPeak;
    metrics.RegionPeakTimeMs(region) = times(index);
    threshold = regionPeak/exp(1);
    recovery = find((1:numel(times))'>index & abs(regional(:,region))<=threshold,1);
    if ~isempty(recovery)
        metrics.RegionPersistenceMs(region) = times(recovery)-times(index);
    end
    metrics.RegionLateResidualDeg(region) = mean(abs(regional(times>=120,region)));
end
metrics.OrthogonalNeighborPersistenceRatio = ...
    metrics.RegionPersistenceMs(3)/metrics.RegionPersistenceMs(2);
end
