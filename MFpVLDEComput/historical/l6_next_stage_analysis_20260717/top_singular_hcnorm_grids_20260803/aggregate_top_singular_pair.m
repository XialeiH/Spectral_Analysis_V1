function aggregate_top_singular_pair()
% Validate and assemble one 51-by-51 top-singular-vector HCnorm scan.

outputRoot = getenv('TSHC_OUTPUT_ROOT');
pairType = string(getenv('TSHC_PAIR'));
combinedFile = fullfile(outputRoot,'all_top_singular_grid_points.tsv');
gridData = readtable(combinedFile,'FileType','text','Delimiter','\t', ...
    'TextType','string');
gridStep = 0.008;
gridMaximum = 0.4;
gridCount = round(gridMaximum/gridStep)+1;
totalCount = gridCount^2;
if height(gridData)~=totalCount || numel(unique(gridData.taskId))~=totalCount
    error('TopSingularHC:Incomplete','Expected %d unique rows, found %d.', ...
        totalCount,height(gridData));
end
gridData = sortrows(gridData,'taskId');
if any(gridData.taskId~=(1:totalCount)')
    error('TopSingularHC:Order','Task ids are incomplete or out of order.');
end

modeHC = nan(gridCount,gridCount);
residual = nan(gridCount,gridCount);
accepted = false(gridCount,gridCount);
gap = nan(gridCount,gridCount);
overlap = nan(gridCount,gridCount);
for row=1:height(gridData)
    ix=gridData.betaXIndex(row);
    iy=gridData.betaYIndex(row);
    modeHC(iy,ix)=gridData.modeHCnorm(row);
    residual(iy,ix)=gridData.fixedPointResidual(row);
    accepted(iy,ix)=logical(gridData.fixedPointAccepted(row)) && ...
        gridData.fixedPointResidual(row)<=1e-4 && ...
        isfinite(gridData.modeHCnorm(row));
    gap(iy,ix)=gridData.singularGapRelative(row);
    overlap(iy,ix)=gridData.modeOverlap(row);
end
validModeHC=modeHC;
validModeHC(~accepted)=NaN;
metrics = table( ...
    ["grid_condition_count";"accepted_condition_count"; ...
    "unaccepted_condition_count";"maximum_accepted_residual"; ...
    "minimum_mode_HCnorm";"maximum_mode_HCnorm"; ...
    "minimum_singular_gap_relative";"minimum_mode_overlap"], ...
    [totalCount;nnz(accepted);nnz(~accepted); ...
    max(residual(accepted),[],'omitnan'); ...
    min(validModeHC,[],'all','omitnan'); ...
    max(validModeHC,[],'all','omitnan'); ...
    min(gap(accepted),[],'omitnan');min(overlap(accepted),[],'omitnan')], ...
    'VariableNames',{'metric','value'});
writetable(metrics,fullfile(outputRoot,'top_singular_grid_metrics.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputRoot,'top_singular_grid_results.mat'),'gridData', ...
    'pairType','modeHC','validModeHC','residual','accepted','gap','overlap', ...
    'metrics','-v7.3');
disp(metrics);
end
