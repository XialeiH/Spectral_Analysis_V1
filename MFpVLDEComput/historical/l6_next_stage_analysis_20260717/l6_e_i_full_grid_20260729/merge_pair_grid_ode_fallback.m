function merge_pair_grid_ode_fallback()
% Merge successful ODE fallbacks into the iteration-only signed grid.

sourceFile = getenv('PAIR_GRID_SOURCE_TABLE');
fallbackFile = getenv('PAIR_GRID_ODE_COMBINED');
outputRoot = getenv('PAIR_GRID_OUTPUT_ROOT');
acceptanceTolerance = str2double(getenv('PAIR_GRID_ODE_ACCEPT_TOL'));
if ~isfinite(acceptanceTolerance) || acceptanceTolerance<=0
    acceptanceTolerance = 1e-9;
end
gridData = readtable(sourceFile,'FileType','text','Delimiter','\t', ...
    'TextType','string');
fallback = readtable(fallbackFile,'FileType','text','Delimiter','\t', ...
    'TextType','string');

gridData.solverMethod = repmat("relaxed_iteration",height(gridData),1);
gridData.odeIntrinsicTime = nan(height(gridData),1);
gridData.odeAcceptedSteps = nan(height(gridData),1);
gridData.odeTermination = repmat("not_attempted",height(gridData),1);
gridData.odeAccepted = false(height(gridData),1);
for index = 1:height(fallback)
    row = find(gridData.taskId==fallback.taskId(index),1);
    if isempty(row) || gridData.converged(row)
        error('PairGridODE:Merge','Invalid fallback task id %d.', ...
            fallback.taskId(index));
    end
    gridData.odeIntrinsicTime(row) = fallback.odeIntrinsicTime(index);
    gridData.odeAcceptedSteps(row) = fallback.odeAcceptedSteps(index);
    gridData.odeTermination(row) = fallback.termination(index);
    accepted = isfinite(fallback.fixedPointResidual(index)) && ...
        fallback.fixedPointResidual(index)<=acceptanceTolerance;
    gridData.odeAccepted(row) = accepted;
    if accepted
        gridData.converged(row) = true;
        gridData.fixedPointResidual(row) = fallback.fixedPointResidual(index);
        gridData.HCnorm(row) = fallback.HCnorm(index);
        gridData.minimumRateHz(row) = fallback.minimumRateHz(index);
        gridData.maximumRateHz(row) = fallback.maximumRateHz(index);
        gridData.solverMethod(row) = "ode45_fallback";
    else
        gridData.solverMethod(row) = "iteration_plus_ode45_failed";
    end
end
writetable(gridData,fullfile(outputRoot,'all_pair_grid_points.tsv'), ...
    'FileType','text','Delimiter','\t');

accepted = isfinite(fallback.fixedPointResidual) & ...
    fallback.fixedPointResidual<=acceptanceTolerance;
solverSummary = table(height(fallback),nnz(accepted),nnz(~accepted), ...
    acceptanceTolerance,max(fallback.odeIntrinsicTime), ...
    max(fallback.odeAcceptedSteps), ...
    'VariableNames',{'attempted','accepted','notAccepted', ...
    'acceptanceTolerance','maximumIntrinsicTime','maximumAcceptedSteps'});
writetable(solverSummary,fullfile(outputRoot,'ode_fallback_metrics.tsv'), ...
    'FileType','text','Delimiter','\t');
aggregate_pair_grid
end
