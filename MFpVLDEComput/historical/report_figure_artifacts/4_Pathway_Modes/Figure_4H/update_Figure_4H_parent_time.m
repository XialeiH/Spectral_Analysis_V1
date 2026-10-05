function update_Figure_4H_parent_time()
% Correct the parent-CG epoch-to-millisecond conversion without rerunning CG.

artifactRoot = fileparts(mfilename('fullpath'));
dataRoot = fullfile(artifactRoot,'data');
analysisFile = fullfile(dataRoot,'Figure_4H_analysis.mat');
parentFile = fullfile(dataRoot,'Figure_4H_parent_CG.mat');

data = load(analysisFile);
parent = load(parentFile,'parentCGStates','parentCGTimesMs', ...
    'parentCGFixedVector');

relaxationP = 0.33;
data.parentCGIterationMs = relaxationP*data.tauMs;
data.parentCGTimesMs = double(parent.parentCGTimesMs(:).')* ...
    data.parentCGIterationMs;
timeMask = data.parentCGTimesMs<=180;
data.parentCGTimesMs = data.parentCGTimesMs(timeMask);
parentDelta = double(parent.parentCGStates(:,timeMask))- ...
    double(parent.parentCGFixedVector(:));

stateRegions = false(numel(data.qTheta),3);
for region = 1:3
    stateRegions(:,region) = repmat(data.spatialRegion(:)==region,3,1);
end
data.parentCGDecoder = (data.qTheta'*parentDelta/data.decoderDenominator).';
data.parentCGRegional = zeros(size(parentDelta,2),3);
for region = 1:3
    maskedQ = data.qTheta.*stateRegions(:,region);
    data.parentCGRegional(:,region) = ...
        (maskedQ'*parentDelta/data.decoderDenominator).';
end
data.parentCGMetrics = local_metrics( ...
    data.parentCGTimesMs,data.parentCGDecoder,data.parentCGRegional);

save(analysisFile,'-struct','data','-v7.3');
fprintf('Parent-CG epoch duration: %.9f ms. Peak at %.6f ms.\n', ...
    data.parentCGIterationMs,data.parentCGMetrics.PeakTimeMs);
fprintf('Parent-CG regional persistence (ms): %s.\n', ...
    mat2str(data.parentCGMetrics.RegionPersistenceMs,6));
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
    recovery = find((1:numel(times))'>index & ...
        abs(regional(:,region))<=threshold,1);
    if ~isempty(recovery)
        metrics.RegionPersistenceMs(region) = times(recovery)-times(index);
    end
    metrics.RegionLateResidualDeg(region) = ...
        mean(abs(regional(times>=120,region)));
end
metrics.OrthogonalNeighborPersistenceRatio = ...
    metrics.RegionPersistenceMs(3)/metrics.RegionPersistenceMs(2);
end
