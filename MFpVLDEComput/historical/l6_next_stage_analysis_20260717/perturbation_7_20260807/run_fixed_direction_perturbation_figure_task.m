function output = run_fixed_direction_perturbation_figure_task(taskId,setupFile,outputRoot)
% Generate one full-format perturbation figure near the HC basin boundary.

arguments
    taskId (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    outputRoot (1,:) char
end

hcValues = [6.00 6.01 6.02 6.03 6.10 6.50];
if taskId>numel(hcValues)
    error('Perturbation7:TaskId','Task %d exceeds %d cases.', ...
        taskId,numel(hcValues));
end
figureNumber = sprintf('7.%d',taskId+5);
requestedHC = hcValues(taskId);
output = run_perturbation_figure_task(4,setupFile,outputRoot, ...
    requestedHC,figureNumber,7006);

tag = strrep(sprintf('%.2f',requestedHC),'.','p');
stem = sprintf(['%s_Perturbation_Transient_and_Return_' ...
    'Baseline_L6_FixedDirection_HC%s'],figureNumber,tag);
replot_perturbation_e_only(fullfile(outputRoot,[stem '_data.mat']), ...
    fullfile(outputRoot,[stem '.pdf']));
end
