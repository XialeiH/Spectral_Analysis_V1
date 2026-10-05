function run_figure6_0_grid_task()
% Compute and save all Figure 6.0 quantities for one contrast/grid point.

taskId = local_task_id();
contrasts = [19 42 66 100];
betaGrid = 0:0.008:0.4;
gridCount = numel(betaGrid);
pointsPerContrast = gridCount^2;
totalCount = numel(contrasts)*pointsPerContrast;
if ~isfinite(taskId) || taskId~=round(taskId) || taskId<1 || taskId>totalCount
    error('Figure60:TaskId','Task id must be an integer in 1:%d.',totalCount);
end

setupRoot = getenv('FIG60_SETUP_ROOT');
outputRoot = getenv('FIG60_OUTPUT_ROOT');
fixedTolerance = str2double(getenv('FIG60_FIXED_TOLERANCE'));
if ~isfinite(fixedTolerance); fixedTolerance=1e-4; end
if isempty(setupRoot) || isempty(outputRoot)
    error('Figure60:Environment','FIG60 setup and output roots are required.');
end

contrastIndex = floor((taskId-1)/pointsPerContrast)+1;
pointId = mod(taskId-1,pointsPerContrast)+1;
beta6Index = floor((pointId-1)/gridCount)+1;
betaIIndex = mod(pointId-1,gridCount)+1;
contrast = contrasts(contrastIndex);
beta6 = betaGrid(beta6Index);
betaI = betaGrid(betaIIndex);
gain6 = 1+beta6;
gainI = 1+betaI;

setupFile = fullfile(setupRoot,sprintf( ...
    'figure6_0_baseline_contrast%d.mat',contrast));
loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
baseline = setup.BaselineState(:);
n = numel(baseline)/3;

phi = @(state)l6ns_phi(state,1-gain6,context,[1 1],gainI,'true');
fixed = real_tuning_fixed_point_tolerance( ...
    phi,baseline,context.RelaxationP,fixedTolerance);
isBaselinePoint = beta6==0 && betaI==0;
if isBaselinePoint
    fixed.State = baseline;
    fixed.Converged = true;
    fixed.Iterations = 0;
    fixed.Residual = max(abs(phi(baseline)-baseline));
end

result = struct();
result.TaskId = taskId;
result.ContrastIndex = contrastIndex;
result.Contrast = contrast;
result.PointId = pointId;
result.Beta6Index = beta6Index;
result.BetaIIndex = betaIIndex;
result.Beta6 = beta6;
result.BetaI = betaI;
result.FixedPointAccepted = fixed.Converged;
result.FixedPointResidual = fixed.Residual;
result.FixedPointIterations = fixed.Iterations;
result.FixedPointHC = NaN;
result.SingularModeHC = NaN;
result.SingularValue = NaN;
result.SingularGapRelative = NaN;
result.SingularResidual = NaN;
result.EigenclusterHC = NaN;
result.ClusterCount = NaN;
result.ClusterRank = NaN;
result.ClusterMaximumRealEigenvalue = NaN;
result.ClusterMaximumResidual = NaN;
result.TuningWassersteinCircular = nan(1,2);
result.TuningWassersteinLinear = nan(1,2);
result.CanonicalFixedPointAccepted = false(1,numel(setup.CanonicalAngles));
result.CanonicalFixedPointResiduals = nan(1,numel(setup.CanonicalAngles));
result.FixedPointState = single(nan(size(baseline)));
result.LeadingRightSingularMode = single(nan(size(baseline)));
result.EigenclusterEnvelope = single(nan(size(baseline)));
result.EigenclusterEigenvalues = complex(single([]));
result.TuningCurves = single(nan(numel(setup.FullAngles),2));

if fixed.Converged
    result.FixedPointState = single(fixed.State);
    result.FixedPointHC = local_hc(fixed.State,baseline,n);
    J = real_tuning_true_jacobian(fixed.State,context,gain6,gainI);

    [rightMode,singularValues,singularResidual] = ...
        local_leading_right_singular_mode(J,setup.BaselineRightMode);
    result.LeadingRightSingularMode = single(rightMode);
    result.SingularValue = singularValues(1);
    result.SingularGapRelative = ...
        (singularValues(1)-singularValues(2))/max(singularValues(1),eps);
    result.SingularResidual = singularResidual;
    result.SingularModeHC = local_hc( ...
        rightMode,setup.BaselineRightMode,n);

    [clusterEnvelope,clusterValues,clusterRank,clusterResidual] = ...
        local_leading_cluster(J,setup.ClusterTolerance);
    result.EigenclusterEnvelope = single(clusterEnvelope);
    result.EigenclusterEigenvalues = single(clusterValues);
    result.ClusterCount = numel(clusterValues);
    result.ClusterRank = clusterRank;
    result.ClusterMaximumRealEigenvalue = max(real(clusterValues));
    result.ClusterMaximumResidual = clusterResidual;
    result.EigenclusterHC = local_hc( ...
        clusterEnvelope,setup.BaselineClusterEnvelope,n);

    canonicalStates = zeros(numel(baseline),numel(setup.CanonicalAngles));
    if isBaselinePoint
        canonicalStates = setup.BaselineCanonicalStates;
        result.CanonicalFixedPointAccepted(:) = true;
        result.CanonicalFixedPointResiduals(:) = 0;
    else
        for angleIndex = 1:numel(setup.CanonicalAngles)
            angleContext = context;
            angleContext.OrientationUse = setup.CanonicalAngles(angleIndex)* ...
                ones(size(context.OrientationUse));
            anglePhi = @(state)l6ns_phi( ...
                state,1-gain6,angleContext,[1 1],gainI,'true');
            angleFixed = real_tuning_fixed_point_tolerance( ...
                anglePhi,setup.BaselineCanonicalStates(:,angleIndex), ...
                angleContext.RelaxationP,fixedTolerance);
            result.CanonicalFixedPointAccepted(angleIndex) = angleFixed.Converged;
            result.CanonicalFixedPointResiduals(angleIndex) = angleFixed.Residual;
            if ~angleFixed.Converged; break; end
            canonicalStates(:,angleIndex) = angleFixed.State;
        end
    end
    if all(result.CanonicalFixedPointAccepted)
        curves = figure6_0_tuning_curves( ...
            canonicalStates,setup.CanonicalAngles,setup.CWeight);
        result.TuningCurves = single(curves);
        for pixelIndex = 1:2
            [result.TuningWassersteinCircular(pixelIndex), ...
                result.TuningWassersteinLinear(pixelIndex)] = ...
                figure6_0_circular_wasserstein( ...
                setup.BaselineTuningCurves(:,pixelIndex), ...
                curves(:,pixelIndex),setup.FullAngles,0.1);
        end
    end
    if isBaselinePoint
        result.FixedPointHC = 0;
        result.SingularModeHC = 0;
        result.EigenclusterHC = 0;
        result.TuningWassersteinCircular(:) = 0;
        result.TuningWassersteinLinear(:) = 0;
    end
end

pointsRoot = fullfile(outputRoot,sprintf('contrast_%03d',contrast),'points');
if ~exist(pointsRoot,'dir'); mkdir(pointsRoot); end
outputFile = fullfile(pointsRoot,sprintf('point_%04d.mat',pointId));
temporaryFile = [outputFile '.tmp.mat'];
save(temporaryFile,'result','-v7');
movefile(temporaryFile,outputFile,'f');
fprintf(['Figure 6.0 %d/%d: contrast %d point %d beta6 %.3f betaI %.3f ' ...
    'accepted %d residual %.3g HC [%.4g %.4g %.4g] W1 [%.4g %.4g].\n'], ...
    taskId,totalCount,contrast,pointId,beta6,betaI,fixed.Converged, ...
    fixed.Residual,result.FixedPointHC,result.SingularModeHC, ...
    result.EigenclusterHC,result.TuningWassersteinCircular);
end

function taskId = local_task_id()
taskId = str2double(getenv('FIG60_TASK_ID'));
arrayId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isfinite(arrayId)
    offset = str2double(getenv('FIG60_TASK_OFFSET'));
    if ~isfinite(offset); offset=0; end
    taskId = arrayId+offset;
end
end

function value = local_hc(state,baseline,n)
value = HC_norm_diff(state(1:n),state(n+(1:n)),state(2*n+(1:n)), ...
    baseline(1:n),baseline(n+(1:n)),baseline(2*n+(1:n)), ...
    0.3077,0.8,0.2);
end

function [rightMode,values,residual] = ...
        local_leading_right_singular_mode(J,baselineMode)
options = struct('tol',1e-7,'maxit',2000,'p',80,'disp',0, ...
    'v0',baselineMode(:));
[leftVectors,S,rightVectors,flag] = svds(J,2,'largest',options);
if flag~=0; error('Figure60:Svds','svds returned flag %d.',flag); end
values = real(diag(S));
[values,order] = sort(values,'descend');
leftVectors = leftVectors(:,order);
rightVectors = rightVectors(:,order);
[~,pivot] = max(abs(rightVectors(:,1)));
phase = exp(-1i*angle(rightVectors(pivot,1)));
rightMode = real(rightVectors(:,1)*phase);
leftMode = real(leftVectors(:,1)*phase);
rightMode = rightMode/max(norm(rightMode),eps);
leftMode = leftMode/max(norm(leftMode),eps);
if dot(rightMode,baselineMode)<0
    rightMode=-rightMode;
    leftMode=-leftMode;
end
residual = max(norm(J*rightMode-values(1)*leftMode), ...
    norm(J'*leftMode-values(1)*rightMode))/max(values(1),eps);
end

function [envelope,selectedValues,clusterRank,maxResidual] = ...
        local_leading_cluster(J,tolerance)
dimension = size(J,1);
requested = 48;
maximum = 192;
while true
    options = struct('tol',1e-9,'maxit',4000, ...
        'p',min(dimension,max(96,2*requested+16)), ...
        'disp',0,'isreal',true);
    [vectors,D,flag] = eigs(J,requested,'largestreal',options);
    if flag~=0; error('Figure60:Eigs','eigs returned flag %d.',flag); end
    values = diag(D);
    [~,order] = sortrows([real(values),imag(values)],[-1 -2]);
    values = values(order);
    vectors = vectors(:,order);
    topReal = real(values(1));
    padding = 100*eps(max(abs(topReal),1));
    selected = find(real(values)>=topReal-(tolerance+padding));
    boundaryIsClear = real(values(end))<topReal-(tolerance+padding);
    if boundaryIsClear || requested>=maximum; break; end
    requested = min(maximum,2*requested);
end
if ~boundaryIsClear
    error('Figure60:ClusterBoundary', ...
        'Top cluster is not closed after %d eigenvalues.',requested);
end
selectedValues = values(selected);
selectedVectors = vectors(:,selected);
selectedVectors = selectedVectors./max(vecnorm(selectedVectors),eps);
[basis,R] = qr(selectedVectors,0);
rankTolerance = max(size(R))*eps(max(norm(R,2),1));
clusterRank = sum(abs(diag(R))>rankTolerance);
basis = basis(:,1:clusterRank);
envelope = sqrt(sum(abs(basis).^2,2));
residuals = vecnorm(J*selectedVectors- ...
    selectedVectors.*reshape(selectedValues,1,[])) ./ ...
    max((norm(J,'fro')+abs(selectedValues.')).*vecnorm(selectedVectors),eps);
maxResidual = max(residuals);
end
