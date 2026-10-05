function aggregate_dense_finite_time_alignment(inputRoot,trajectoryRoot,outputRoot)
% Combine dense finite-time array outputs into one curve per figure.

sampleCount = 251;
figureCount = 12;
if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
for figureIndex = 1:figureCount
    trajectory = load(fullfile(trajectoryRoot, ...
        sprintf('figure_%02d_trajectory.mat',figureIndex)),'sampleTimesMs');
    figureSampleCount = numel(trajectory.sampleTimesMs);
    timesMs = zeros(figureSampleCount,1);
    alignment = zeros(figureSampleCount,1);
    peakHorizonMs = zeros(figureSampleCount,1);
    peakGain = zeros(figureSampleCount,1);
    extendedTo100Ms = false(figureSampleCount,1);
    for sampleIndex = 1:figureSampleCount
        if sampleIndex<=sampleCount
            taskId = (figureIndex-1)*sampleCount+sampleIndex;
            inputFile = fullfile(inputRoot,sprintf('task_%04d.mat',taskId));
        else
            inputFile = fullfile(inputRoot,sprintf( ...
                'figure_%02d_sample_%03d.mat',figureIndex,sampleIndex));
        end
        loaded = load(inputFile,'result');
        current = loaded.result;
        timesMs(sampleIndex) = current.TimeMs;
        alignment(sampleIndex) = current.Alignment;
        peakHorizonMs(sampleIndex) = current.PeakHorizonMs;
        peakGain(sampleIndex) = current.PeakGain;
        extendedTo100Ms(sampleIndex) = current.ExtendedTo100Ms;
    end
    curve = table(timesMs,alignment,peakHorizonMs,peakGain,extendedTo100Ms);
    save(fullfile(outputRoot,sprintf('figure_%02d_curve.mat',figureIndex)),'curve');
    writetable(curve,fullfile(outputRoot,sprintf('figure_%02d_curve.tsv',figureIndex)), ...
        'FileType','text','Delimiter','\t');
end
end
