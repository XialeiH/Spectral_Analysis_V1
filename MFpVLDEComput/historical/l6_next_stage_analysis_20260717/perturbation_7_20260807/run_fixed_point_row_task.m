function run_fixed_point_row_task(taskId,setupFile,dataRoot,outputRoot)
% Add a fixed-point row to HC 6.00, 6.01, or 6.02.

arguments
    taskId (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    dataRoot (1,:) char
    outputRoot (1,:) char
end

hcTags = {'6p00','6p01','6p02'};
figureNumbers = {'7.6','7.7','7.8'};
if taskId>numel(hcTags)
    error('Perturbation7:TaskId','Task %d exceeds %d cases.',taskId,numel(hcTags));
end
stem = sprintf(['%s_Perturbation_Transient_and_Return_' ...
    'Baseline_L6_FixedDirection_HC%s'],figureNumbers{taskId},hcTags{taskId});
if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
replot_with_fixed_point_row(fullfile(dataRoot,[stem '_data.mat']),setupFile, ...
    fullfile(outputRoot,[stem '.pdf']), ...
    fullfile(outputRoot,[stem '_fixedpointrow_data.mat']));
end
