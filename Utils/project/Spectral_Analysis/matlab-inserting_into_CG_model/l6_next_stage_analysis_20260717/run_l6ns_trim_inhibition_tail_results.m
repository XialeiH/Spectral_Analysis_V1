% Truncate repeated inhibition laps at their closest return to the persistent branch.

sourceRoot = getenv('L6NS_GLOBAL_BIFURCATION_ARCLENGTH_ROOT');
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_FINAL_BRANCH_ROOT');
if isempty(sourceRoot) || isempty(outputDir)
    error('Source arclength and final branch roots are required.');
end
for branchSign = [-1 1]
    stem = sprintf('inhibition_sign_%+d',branchSign);
    stem = strrep(stem,'+','p');
    stem = strrep(stem,'-','m');
    sourceLoaded = load(fullfile(sourceRoot,[stem '_result.mat']),'result');
    finalLoaded = load(fullfile(outputDir,[stem '_result.mat']),'result');
    sourceDistance = sourceLoaded.result.Table.branchDistance;
    jump = find(sourceDistance(2:end)<0.2*sourceDistance(1:end-1) & ...
        sourceDistance(1:end-1)>1,1,'first');
    if isempty(jump)
        error('No source branch-switch jump found for sign %+d.',branchSign);
    end
    result = finalLoaded.result;
    appendedDistance = result.Table.branchDistance(jump+1:end);
    [closestDistance,relativeIndex] = min(appendedDistance);
    keep = jump+relativeIndex;
    result.Table = result.Table(1:keep,:);
    result.States = result.States(1:keep);
    result.TerminalReason = "closest_return_to_persistent";
    writetable(result.Table,fullfile(outputDir,[stem '_continuation.tsv']), ...
        'FileType','text','Delimiter','\t');
    save(fullfile(outputDir,[stem '_result.mat']),'result','-v7.3');
    fprintf('Inhibition sign %+d trimmed at row %d: w=%.12g distance=%.6g amplitude=%.6g\n', ...
        branchSign,keep,result.Table.freezeWeight(end),closestDistance, ...
        result.Table.normalizedAmplitude(end));
end
