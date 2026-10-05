function aggregate_figure6_0_contrast()
% Assemble one contrast's 51-by-51 Figure 6.0 dataset.

contrast = str2double(getenv('FIG60_CONTRAST'));
outputRoot = getenv('FIG60_OUTPUT_ROOT');
setupRoot = getenv('FIG60_SETUP_ROOT');
if ~ismember(contrast,[19 42 66 100]) || isempty(outputRoot)
    error('Figure60:AggregateEnvironment','Valid contrast and output root required.');
end
loaded = load(fullfile(setupRoot,sprintf( ...
    'figure6_0_baseline_contrast%d.mat',contrast)),'setup');
setup = loaded.setup;

betaGrid = 0:0.008:0.4;
gridCount = numel(betaGrid);
pointCount = gridCount^2;
stateCount = numel(setup.BaselineState);
angleCount = numel(setup.FullAngles);
pointsRoot = fullfile(outputRoot,sprintf('contrast_%03d',contrast),'points');

fixedPointStates = nan(stateCount,pointCount,'single');
leadingRightSingularModes = nan(stateCount,pointCount,'single');
leadingEigenclusterEnvelopes = nan(stateCount,pointCount,'single');
tuningCurves = nan(angleCount,2,pointCount,'single');
fixedPointHC = nan(gridCount);
singularModeHC = nan(gridCount);
eigenclusterHC = nan(gridCount);
tuningWassersteinCircular = nan(gridCount,gridCount,2);
tuningWassersteinLinear = nan(gridCount,gridCount,2);
fixedPointResidual = nan(gridCount);
fixedPointAccepted = false(gridCount);
canonicalAccepted = false(gridCount);
singularValue = nan(gridCount);
singularGapRelative = nan(gridCount);
clusterCount = nan(gridCount);
clusterRank = nan(gridCount);
clusterMaximumRealEigenvalue = nan(gridCount);
clusterEigenvalues = cell(pointCount,1);

summaryRows = cell(pointCount,22);
for pointId = 1:pointCount
    pointFile = fullfile(pointsRoot,sprintf('point_%04d.mat',pointId));
    if ~isfile(pointFile)
        error('Figure60:MissingPoint','Missing %s.',pointFile);
    end
    loadedPoint = load(pointFile,'result');
    result = loadedPoint.result;
    if result.PointId~=pointId || result.Contrast~=contrast
        error('Figure60:PointIdentity','Point %d metadata mismatch.',pointId);
    end
    ix = result.Beta6Index;
    iy = result.BetaIIndex;
    fixedPointStates(:,pointId) = result.FixedPointState;
    leadingRightSingularModes(:,pointId) = result.LeadingRightSingularMode;
    leadingEigenclusterEnvelopes(:,pointId) = result.EigenclusterEnvelope;
    tuningCurves(:,:,pointId) = result.TuningCurves;
    clusterEigenvalues{pointId} = result.EigenclusterEigenvalues;
    fixedPointHC(iy,ix) = result.FixedPointHC;
    singularModeHC(iy,ix) = result.SingularModeHC;
    eigenclusterHC(iy,ix) = result.EigenclusterHC;
    tuningWassersteinCircular(iy,ix,:) = ...
        reshape(result.TuningWassersteinCircular,1,1,2);
    tuningWassersteinLinear(iy,ix,:) = ...
        reshape(result.TuningWassersteinLinear,1,1,2);
    fixedPointResidual(iy,ix) = result.FixedPointResidual;
    fixedPointAccepted(iy,ix) = result.FixedPointAccepted;
    canonicalAccepted(iy,ix) = all(result.CanonicalFixedPointAccepted);
    singularValue(iy,ix) = result.SingularValue;
    singularGapRelative(iy,ix) = result.SingularGapRelative;
    clusterCount(iy,ix) = result.ClusterCount;
    clusterRank(iy,ix) = result.ClusterRank;
    clusterMaximumRealEigenvalue(iy,ix) = ...
        result.ClusterMaximumRealEigenvalue;
    summaryRows(pointId,:) = {result.TaskId,contrast,pointId,ix,iy, ...
        result.Beta6,result.BetaI,result.FixedPointAccepted, ...
        result.FixedPointResidual,result.FixedPointIterations, ...
        all(result.CanonicalFixedPointAccepted), ...
        max(result.CanonicalFixedPointResiduals,[],'omitnan'), ...
        result.FixedPointHC,result.SingularModeHC,result.SingularValue, ...
        result.SingularGapRelative,result.SingularResidual, ...
        result.EigenclusterHC,result.ClusterCount,result.ClusterRank, ...
        result.TuningWassersteinCircular(1), ...
        result.TuningWassersteinCircular(2)};
end

pointSummary = cell2table(summaryRows,'VariableNames',{ ...
    'taskId','contrast','pointId','beta6Index','betaIIndex','beta6','betaI', ...
    'fixedPointAccepted','fixedPointResidual','fixedPointIterations', ...
    'canonicalFixedPointsAccepted','maximumCanonicalResidual', ...
    'fixedPointHC','singularModeHC','singularValue', ...
    'singularGapRelative','singularResidual','eigenclusterHC', ...
    'clusterCount','clusterRank','pixel5x10CircularW1Deg', ...
    'pixel1x10CircularW1Deg'});

metricNames = {'fixedPointHC','singularModeHC','eigenclusterHC', ...
    'pixel5x10CircularW1Deg','pixel1x10CircularW1Deg'};
metricMaps = cat(3,fixedPointHC,singularModeHC,eigenclusterHC, ...
    tuningWassersteinCircular(:,:,1),tuningWassersteinCircular(:,:,2));
regressionSlope = nan(numel(metricNames),1);
regressionR2 = nan(numel(metricNames),1);
for metricIndex = 1:numel(metricNames)
    [regressionSlope(metricIndex),regressionR2(metricIndex)] = ...
        local_minimum_path_regression(betaGrid,metricMaps(:,:,metricIndex));
end
regression = table(string(metricNames(:)),regressionSlope,regressionR2, ...
    'VariableNames',{'metric','throughOriginSlopeBetaIOverBeta6','R2'});

contrastRoot = fullfile(outputRoot,sprintf('contrast_%03d',contrast));
writetable(pointSummary,fullfile(contrastRoot,'point_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(regression,fullfile(contrastRoot,'metric_regression.tsv'), ...
    'FileType','text','Delimiter','\t');
datasetFile = fullfile(contrastRoot,sprintf( ...
    'figure6_0_contrast%d_dataset.mat',contrast));
save(datasetFile,'contrast','betaGrid','setup','pointSummary','regression', ...
    'fixedPointStates','leadingRightSingularModes', ...
    'leadingEigenclusterEnvelopes','clusterEigenvalues','tuningCurves', ...
    'fixedPointHC','singularModeHC','eigenclusterHC', ...
    'tuningWassersteinCircular','tuningWassersteinLinear', ...
    'fixedPointResidual','fixedPointAccepted','canonicalAccepted', ...
    'singularValue','singularGapRelative','clusterCount','clusterRank', ...
    'clusterMaximumRealEigenvalue','-v7.3');
fprintf('Aggregated Figure 6.0 contrast %d: %d points -> %s\n', ...
    contrast,pointCount,datasetFile);
end

function [slope,rSquared] = local_minimum_path_regression(betaGrid,metric)
bestBetaI = nan(size(betaGrid));
for column = 1:numel(betaGrid)
    values = metric(:,column);
    finiteRows = find(isfinite(values));
    if isempty(finiteRows); continue; end
    [~,localIndex] = min(values(finiteRows));
    bestBetaI(column) = betaGrid(finiteRows(localIndex));
end
valid = isfinite(bestBetaI) & betaGrid>0;
x = betaGrid(valid).';
y = bestBetaI(valid).';
if isempty(x); slope=NaN; rSquared=NaN; return; end
slope = (x'*y)/(x'*x);
prediction = slope*x;
rSquared = 1-sum((y-prediction).^2)/max(sum((y-mean(y)).^2),eps);
end
