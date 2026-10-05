function compute_dense_finite_time_alignment_extra_block_task(blockId,setupFile,trajectoryRoot,outputRoot)
% Compute three samples from t=251:500 ms for Figures 7.3-7.12.

arguments
    blockId (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    trajectoryRoot (1,:) char
    outputRoot (1,:) char
end

samplesPerFigure = 250;
figureCount = 10;
samplesPerBlock = 3;
totalSamples = samplesPerFigure*figureCount;
firstExtraId = (blockId-1)*samplesPerBlock+1;
lastExtraId = min(firstExtraId+samplesPerBlock-1,totalSamples);
if firstExtraId>totalSamples
    error('Block ID %d exceeds the required range.',blockId);
end
for extraId = firstExtraId:lastExtraId
    figureIndex = floor((extraId-1)/samplesPerFigure)+3;
    sampleIndex = mod(extraId-1,samplesPerFigure)+252;
    compute_dense_finite_time_alignment_task( ...
        extraId,setupFile,trajectoryRoot,outputRoot,figureIndex,sampleIndex);
end
end
