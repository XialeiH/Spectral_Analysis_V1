function prepare_figure6_0_contrast_setups()
% Build authoritative baseline data for the four Figure 6.0 contrasts.

templateFile = getenv('FIG60_TEMPLATE_SETUP');
geometryRoot = getenv('FIG60_GEOMETRY_ROOT');
atlasRoot = getenv('FIG60_ATLAS_ROOT');
outputRoot = getenv('FIG60_SETUP_ROOT');
if ~isfile(templateFile) || ~isfolder(geometryRoot) || ...
        ~isfolder(atlasRoot) || isempty(outputRoot)
    error('Figure60:PrepareEnvironment','Required setup paths are missing.');
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

template = load(templateFile,'setup');
contrasts = [19 42 66 100];
canonicalAngles = 0:3.75:22.5;
for contrast = contrasts
    geometryFile = fullfile(geometryRoot,sprintf( ...
        'geometry_sections4_5_h96baseline_L6eqW0p00_contr%d_angle_0.00.mat', ...
        contrast));
    atlasFile = fullfile(atlasRoot,sprintf( ...
        'full_eigenspectrum_maps_angle0p00_contrast%d.mat',contrast));
    if ~isfile(geometryFile) || ~isfile(atlasFile)
        error('Figure60:PrepareInput','Missing contrast-%d source data.',contrast);
    end
    geometry = load(geometryFile,'Section4','GeometryMetrics');
    atlas = load(atlasFile,'conditionData');
    maps = geometry.GeometryMetrics.FixedPointMaps;
    baselineState = [maps.S(:);maps.C(:);maps.I(:)];
    context = template.setup.Context;
    context.FixedPoint = baselineState;
    context.ContrastUse = contrast*ones(size(context.ContrastUse));
    context.OrientationUse = zeros(size(context.OrientationUse));
    jBaseline = sparse(geometry.Section4.A);
    if ~isequal(size(jBaseline),[4800 4800]) || numel(baselineState)~=4800
        error('Figure60:PrepareDimension','Unexpected contrast-%d dimensions.',contrast);
    end

    options = struct('tol',1e-9,'maxit',2500,'p',80,'disp',0);
    [~,singularMatrix,rightVectors,flag] = svds(jBaseline,2,'largest',options);
    if flag~=0; error('Figure60:PrepareSvds','svds flag %d.',flag); end
    singularValues = diag(singularMatrix);
    [singularValues,order] = sort(real(singularValues),'descend');
    baselineRightMode = local_real_vector(rightVectors(:,order(1)));

    clusterMaps = atlas.conditionData.Jacobians.J_full.TopClusterEnvelopeMaps;
    baselineClusterEnvelope = [clusterMaps.S(:);clusterMaps.C(:);clusterMaps.I(:)];
    baselineClusterCount = ...
        atlas.conditionData.Jacobians.J_full.TopClusterCount;
    baselineClusterRank = ...
        atlas.conditionData.Jacobians.J_full.TopClusterRank;

    baselineCanonicalStates = zeros(4800,numel(canonicalAngles));
    baselineCanonicalStates(:,1) = baselineState;
    previousState = baselineState;
    for angleIndex = 2:numel(canonicalAngles)
        angle = canonicalAngles(angleIndex);
        angleContext = context;
        angleContext.OrientationUse = angle*ones(size(context.OrientationUse));
        phi = @(state)l6ns_phi(state,0,angleContext,[1 1],1,'true');
        fixed = real_tuning_fixed_point_tolerance( ...
            phi,previousState,angleContext.RelaxationP,1e-6);
        if ~fixed.Converged
            error('Figure60:PrepareFixedPoint', ...
                'Baseline contrast %d angle %.2f residual %.6g.', ...
                contrast,angle,fixed.Residual);
        end
        baselineCanonicalStates(:,angleIndex) = fixed.State;
        previousState = fixed.State;
    end
    baselineTuningCurves = figure6_0_tuning_curves( ...
        baselineCanonicalStates,canonicalAngles,context.CWeight);

    setup = struct();
    setup.Contrast = contrast;
    setup.Context = context;
    setup.BaselineState = baselineState;
    setup.BaselineRightMode = baselineRightMode;
    setup.BaselineSingularValues = singularValues(:);
    setup.BaselineClusterEnvelope = baselineClusterEnvelope;
    setup.BaselineClusterCount = baselineClusterCount;
    setup.BaselineClusterRank = baselineClusterRank;
    setup.BaselineCanonicalStates = baselineCanonicalStates;
    setup.BaselineTuningCurves = baselineTuningCurves;
    setup.CanonicalAngles = canonicalAngles;
    setup.FullAngles = 0:3.75:180;
    setup.ClusterTolerance = 1e-3;
    setup.CWeight = context.CWeight;
    setup.SourceGeometry = geometryFile;
    setup.SourceAtlas = atlasFile;
    outputFile = fullfile(outputRoot,sprintf( ...
        'figure6_0_baseline_contrast%d.mat',contrast));
    save(outputFile,'setup','-v7.3');
    fprintf(['Prepared contrast %d: sigma1 %.8g, cluster count %d, ' ...
        'rank %d -> %s\n'],contrast,singularValues(1), ...
        baselineClusterCount,baselineClusterRank,outputFile);
end
end

function vector = local_real_vector(vector)
[~,pivot] = max(abs(vector));
vector = real(vector*exp(-1i*angle(vector(pivot))));
vector = vector/max(norm(vector),eps);
end
