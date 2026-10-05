function compute_dense_finite_time_alignment_block_task(blockId,setupFile,trajectoryRoot,outputRoot)
% Compute three consecutive dense finite-time samples in one Slurm task.

arguments
    blockId (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    trajectoryRoot (1,:) char
    outputRoot (1,:) char
end

samplesPerBlock = 3;
totalSamples = 12*251;
firstTaskId = (blockId-1)*samplesPerBlock+1;
lastTaskId = min(firstTaskId+samplesPerBlock-1,totalSamples);
if firstTaskId>totalSamples
    error('Block ID %d exceeds the required range.',blockId);
end
for taskId = firstTaskId:lastTaskId
    compute_dense_finite_time_alignment_task( ...
        taskId,setupFile,trajectoryRoot,outputRoot);
end
end
