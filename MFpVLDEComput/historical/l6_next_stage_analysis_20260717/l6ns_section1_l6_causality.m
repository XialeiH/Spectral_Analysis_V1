function result = l6ns_section1_l6_causality(cfg, data)
% Section 1: establish the causal L6 origin of the positive branch.

sectionDir = fullfile(cfg.OutputRoot, 'section1_l6_causality');
if ~exist(sectionDir, 'dir'); mkdir(sectionDir); end
rng(cfg.RandomSeed);

% Direct boundary and generalized-pencil prediction.
directBoundary = l6ns_direct_boundary(data.A, data.B, cfg);
pencilThreshold = l6ns_generalized_threshold( ...
    data.A, data.B, cfg, directBoundary.Alpha);

% Dense, biorthogonally tracked critical branch.
wGrid = cfg.DenseWGrid(:);
count = numel(wGrid);
lambda = nan(count, 1);
kappa = nan(count, 1);
rightResidual = nan(count, 1);
leftResidual = nan(count, 1);
cBase = nan(count, 1);
cL6 = nan(count, 1);
sensitivityAlpha = nan(count, 1);
radialSensitivityAlpha = nan(count, 1);
selectedRight = cell(count, 1);
selectedLeft = cell(count, 1);
matchScore = nan(count, 1);

[~, anchorIndex] = min(abs(wGrid - directBoundary.W));
anchorModes = l6ns_eigenpairs(data.A + (1-wGrid(anchorIndex))*data.B, ...
    cfg.BranchPoolSize, 'largestreal', cfg);
[~, anchorModeIndex] = min(abs(anchorModes.Lambda-cfg.CriticalBoundary));
[lambda,kappa,rightResidual,leftResidual,cBase,cL6,sensitivityAlpha, ...
    radialSensitivityAlpha,selectedRight,selectedLeft,matchScore] = ...
    local_store_mode(anchorIndex,anchorModeIndex,anchorModes,wGrid,data, ...
    lambda,kappa,rightResidual,leftResidual,cBase,cL6,sensitivityAlpha, ...
    radialSensitivityAlpha,selectedRight,selectedLeft,matchScore,1);

previous = struct('Lambda',lambda(anchorIndex),'Right',selectedRight{anchorIndex}, ...
    'Left',selectedLeft{anchorIndex});
for wi = anchorIndex+1:count
    modes = l6ns_eigenpairs(data.A+(1-wGrid(wi))*data.B, ...
        cfg.BranchPoolSize,'largestreal',cfg);
    [index,score] = local_match_mode(previous,modes);
    [lambda,kappa,rightResidual,leftResidual,cBase,cL6,sensitivityAlpha, ...
        radialSensitivityAlpha,selectedRight,selectedLeft,matchScore] = ...
        local_store_mode(wi,index,modes,wGrid,data,lambda,kappa,rightResidual, ...
        leftResidual,cBase,cL6,sensitivityAlpha,radialSensitivityAlpha, ...
        selectedRight,selectedLeft,matchScore,score);
    previous = struct('Lambda',lambda(wi),'Right',selectedRight{wi},'Left',selectedLeft{wi});
end
previous = struct('Lambda',lambda(anchorIndex),'Right',selectedRight{anchorIndex}, ...
    'Left',selectedLeft{anchorIndex});
for wi = anchorIndex-1:-1:1
    modes = l6ns_eigenpairs(data.A+(1-wGrid(wi))*data.B, ...
        cfg.BranchPoolSize,'largestreal',cfg);
    [index,score] = local_match_mode(previous,modes);
    [lambda,kappa,rightResidual,leftResidual,cBase,cL6,sensitivityAlpha, ...
        radialSensitivityAlpha,selectedRight,selectedLeft,matchScore] = ...
        local_store_mode(wi,index,modes,wGrid,data,lambda,kappa,rightResidual, ...
        leftResidual,cBase,cL6,sensitivityAlpha,radialSensitivityAlpha, ...
        selectedRight,selectedLeft,matchScore,score);
    previous = struct('Lambda',lambda(wi),'Right',selectedRight{wi},'Left',selectedLeft{wi});
end

finiteDifferenceW = gradient(lambda, wGrid);
analyticDerivativeW = -sensitivityAlpha;
derivativeRelativeError = abs(finiteDifferenceW - analyticDerivativeW) ./ ...
    max(abs(analyticDerivativeW), 1e-12);

branchTable = table(wGrid, 1-wGrid, real(lambda), imag(lambda), ...
    real(cBase), imag(cBase), real(cL6), imag(cL6), ...
    real(sensitivityAlpha), imag(sensitivityAlpha), ...
    real(analyticDerivativeW), imag(analyticDerivativeW), ...
    real(finiteDifferenceW), imag(finiteDifferenceW), derivativeRelativeError, ...
    radialSensitivityAlpha, kappa, rightResidual, leftResidual, matchScore, ...
    'VariableNames', {'w','alpha','lambdaReal','lambdaImag', ...
    'cBaseReal','cBaseImag','cL6Real','cL6Imag', ...
    'dLambdaDAlphaReal','dLambdaDAlphaImag','dLambdaDwReal','dLambdaDwImag', ...
    'finiteDifferenceDwReal','finiteDifferenceDwImag','derivativeRelativeError', ...
    'radialDAbsLambdaDAlpha','conditionNumber','rightResidual','leftResidual','matchScore'});
writetable(branchTable, fullfile(sectionDir, 'critical_branch.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

% Negative-w asymptotic scaling and positive-edge subspace inheritance.
bSubspace = l6ns_invariant_subspace(data.B,cfg.PositiveSubspaceSize,'largestreal',cfg);
bModes = bSubspace.Modes;
qB = bSubspace.Basis;
bPositiveEdge = max(real(bModes.Lambda));
asymW = cfg.AsymptoticWGrid(:);
asymRows = nan(numel(asymW), 7);
for wi = 1:numel(asymW)
    w = asymW(wi);
    alpha = 1 - w;
    jSubspace = l6ns_invariant_subspace(data.A+alpha*data.B, ...
        cfg.PositiveSubspaceSize,'largestreal',cfg);
    modes = jSubspace.Modes;
    qJ = jSubspace.Basis;
    cosines = svd(qB' * qJ);
    asymRows(wi, :) = [w, alpha, max(real(modes.Lambda)), ...
        max(real(modes.Lambda))/alpha, median(cosines), min(cosines), max(cosines)];
end
asymptoticTable = array2table(asymRows, 'VariableNames', ...
    {'w','alpha','leadingReal','leadingRealOverAlpha', ...
    'medianPrincipalCosine','minimumPrincipalCosine','maximumPrincipalCosine'});
writetable(asymptoticTable, fullfile(sectionDir, 'negative_w_asymptotic.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

% Critical-mode inheritance and full vector-action decomposition.
criticalModes = l6ns_eigenpairs(data.A+directBoundary.Alpha*data.B, ...
    cfg.BranchPoolSize,'largestreal',cfg);
[~,criticalIndex] = min(abs(criticalModes.Lambda-cfg.CriticalBoundary));
rCritical = criticalModes.Right(:,criticalIndex);
lCritical = criticalModes.Left(:,criticalIndex);
alphaCritical = directBoundary.Alpha;
projectionPositiveB = norm(qB' * rCritical)^2 / max(norm(rCritical)^2, eps);

[~, singularValuesB, activeRightB] = svds(data.B, cfg.ActiveSingularSize, 'largest');
projectionLeadingSingularB = norm(activeRightB' * rCritical)^2 / ...
    max(norm(rCritical)^2, eps);
[rowSpaceCoefficient,rowSpaceFlag,rowSpaceRelativeResidual,rowSpaceIterations] = ...
    lsqr(data.B',rCritical,1e-10,1200);
rowSpaceProjection = data.B' * rowSpaceCoefficient;
projectionActiveRowSpaceB = norm(rowSpaceProjection)^2 / ...
    max(norm(rCritical)^2,eps);
muB = (rCritical' * data.B * rCritical) / (rCritical' * rCritical);
etaB = norm(data.B * rCritical - muB * rCritical) / max(norm(data.B * rCritical), eps);

actionA = data.A * rCritical;
actionB = alphaCritical * data.B * rCritical;
rhoA = rCritical' * actionA;
rhoB = rCritical' * actionB;
uA = actionA - rhoA * rCritical;
uB = actionB - rhoB * rCritical;
transverseCancellation = norm(uA + uB) / max(norm(uA) + norm(uB), eps);

% First-order block attribution and compensated derivative ablations.
[blockTable, ablationTable] = local_block_analysis( ...
    data.A, data.B, lCritical, rCritical, cfg, directBoundary);
writetable(blockTable, fullfile(sectionDir, 'critical_block_sensitivity.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
writetable(ablationTable, fullfile(sectionDir, 'critical_block_ablation.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

% Figures.
local_plot_branch(branchTable, directBoundary, pencilThreshold, sectionDir);
local_plot_asymptotic(asymptoticTable, bPositiveEdge, sectionDir);
local_plot_actions(rCritical, actionA, actionB, uA, uB, data, sectionDir);
local_plot_ablation(ablationTable, sectionDir);

summaryTable = table(data.Angle, data.Contrast, directBoundary.W, directBoundary.Alpha, ...
    pencilThreshold.W, pencilThreshold.Alpha, pencilThreshold.Residual, ...
    abs(directBoundary.W - pencilThreshold.W), ...
    abs(rCritical' * pencilThreshold.RightVector), bPositiveEdge, ...
    projectionPositiveB, projectionLeadingSingularB, projectionActiveRowSpaceB, ...
    rowSpaceFlag, rowSpaceRelativeResidual, rowSpaceIterations, etaB, ...
    real(rhoA), imag(rhoA), real(rhoB), imag(rhoB), ...
    norm(uA), norm(uB), transverseCancellation, ...
    median(derivativeRelativeError(2:end-1), 'omitnan'), ...
    max(kappa), max(rightResidual), max(leftResidual), ...
    'VariableNames', {'angle','contrast','directWc','directAlphaC', ...
    'pencilWc','pencilAlphaC','pencilResidual','thresholdAbsoluteDifference', ...
    'criticalVectorPencilOverlap','positiveEdgeB','projectionPositiveB', ...
    'projectionLeadingRightSingularB','projectionActiveRowSpaceB', ...
    'activeRowSpaceLsqrFlag','activeRowSpaceLsqrRelativeResidual', ...
    'activeRowSpaceLsqrIterations','approximateBEigenResidual', ...
    'rhoAReal','rhoAImag','rhoBReal','rhoBImag','uANorm','uBNorm', ...
    'transverseCancellation','medianDerivativeRelativeError', ...
    'maximumConditionNumber','maximumRightResidual','maximumLeftResidual'});
writetable(summaryTable, fullfile(sectionDir, 'section1_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

result = struct();
result.Summary = summaryTable;
result.Branch = branchTable;
result.Asymptotic = asymptoticTable;
result.BlockSensitivity = blockTable;
result.BlockAblation = ablationTable;
result.CriticalRight = rCritical;
result.CriticalLeft = lCritical;
result.ActiveSingularValuesB = diag(singularValuesB);
save(fullfile(sectionDir, 'section1_result.mat'), 'result', '-v7.3');
end

function [index, score] = local_match_mode(previous, modes)
candidateCount = numel(modes.Lambda);
scores = zeros(candidateCount, 1);
scale = max(0.03, 0.12 * max(1, abs(previous.Lambda)));
for j = 1:candidateCount
    overlapForward = abs(previous.Left' * modes.Right(:,j)) / ...
        max(norm(previous.Left) * norm(modes.Right(:,j)), eps);
    overlapBackward = abs(modes.Left(:,j)' * previous.Right) / ...
        max(norm(modes.Left(:,j)) * norm(previous.Right), eps);
    proximity = exp(-abs(modes.Lambda(j) - previous.Lambda) / scale);
    scores(j) = 0.4 * overlapForward + 0.4 * overlapBackward + 0.2 * proximity;
end
[score, index] = max(scores);
end

function [sensitivityTable, ablationTable] = local_block_analysis(a, b, l, r, cfg, referenceBoundary)
n = size(a,1) / 3;
allRows = 1:3*n;
allCols = 1:3*n;
definitions = { ...
    'target_S', 1:n, allCols; ...
    'target_C', n+(1:n), allCols; ...
    'target_I', 2*n+(1:n), allCols; ...
    'source_S', allRows, 1:n; ...
    'source_C', allRows, n+(1:n)};

names = strings(size(definitions,1)+1,1);
sensitivities = nan(size(names));
wc = nan(size(names));
alphaC = nan(size(names));
thresholdShift = nan(size(names));
boundaryStatus = strings(size(names));
leadingAtReference = nan(size(names));
referenceModeOverlap = nan(size(names));
boundaryModeOverlap = nan(size(names));
names(1) = "none";
sensitivities(1) = l' * b * r;
wc(1) = referenceBoundary.W;
alphaC(1) = referenceBoundary.Alpha;
thresholdShift(1) = 0;
boundaryStatus(1) = "reference";
leadingAtReference(1) = cfg.CriticalBoundary;
referenceModeOverlap(1) = 1;
boundaryModeOverlap(1) = 1;

for i = 1:size(definitions,1)
    name = string(definitions{i,1});
    rows = definitions{i,2};
    cols = definitions{i,3};
    block = sparse(size(b,1), size(b,2));
    block(rows, cols) = b(rows, cols);
    names(i+1) = name;
    sensitivities(i+1) = l' * block * r;

    bAblated = b - block;
    referenceModes=l6ns_eigenpairs(a+referenceBoundary.Alpha*bAblated, ...
        cfg.BranchPoolSize,'largestreal',cfg);
    [leadingAtReference(i+1),referenceModeIndex]=max(real(referenceModes.Lambda));
    referenceModeOverlap(i+1)=abs(referenceModes.Right(:,referenceModeIndex)'*r)/ ...
        max(norm(referenceModes.Right(:,referenceModeIndex))*norm(r),eps);
    try
        boundary = l6ns_direct_boundary(a, bAblated, cfg);
        wc(i+1) = boundary.W;
        alphaC(i+1) = boundary.Alpha;
        thresholdShift(i+1) = boundary.W - referenceBoundary.W;
        boundaryStatus(i+1) = "crossing_found";
        boundaryModes=l6ns_eigenpairs(a+boundary.Alpha*bAblated, ...
            cfg.BranchPoolSize,'largestreal',cfg);
        [~,boundaryModeIndex]=min(abs(boundaryModes.Lambda-cfg.CriticalBoundary));
        boundaryModeOverlap(i+1)=abs(boundaryModes.Right(:,boundaryModeIndex)'*r)/ ...
            max(norm(boundaryModes.Right(:,boundaryModeIndex))*norm(r),eps);
    catch exception
        warning('No compensated boundary for %s: %s',name,exception.message);
        wc(i+1) = NaN;
        alphaC(i+1) = NaN;
        thresholdShift(i+1) = NaN;
        boundaryStatus(i+1) = "no_crossing_in_search_range";
    end
end

sensitivityTable = table(names, real(sensitivities), imag(sensitivities), ...
    'VariableNames', {'block','sensitivityReal','sensitivityImag'});
ablationTable = table(names,wc,alphaC,thresholdShift,leadingAtReference, ...
    referenceModeOverlap,boundaryModeOverlap,boundaryStatus, ...
    'VariableNames', {'removedBlock','wCritical','alphaCritical','wCriticalShift', ...
    'leadingRealAtReferenceWc','referenceModeOverlap','boundaryModeOverlap','status'});
end

function [lambda,kappa,rightResidual,leftResidual,cBase,cL6,sensitivityAlpha, ...
    radialSensitivityAlpha,selectedRight,selectedLeft,matchScore] = ...
    local_store_mode(wi,index,modes,wGrid,data,lambda,kappa,rightResidual, ...
    leftResidual,cBase,cL6,sensitivityAlpha,radialSensitivityAlpha, ...
    selectedRight,selectedLeft,matchScore,score)
w = wGrid(wi);
alpha = 1-w;
r = modes.Right(:,index);
l = modes.Left(:,index);
lambda(wi) = modes.Lambda(index);
kappa(wi) = modes.ConditionNumber(index);
rightResidual(wi) = modes.RightResidual(index);
leftResidual(wi) = modes.LeftResidual(index);
cBase(wi) = l'*data.A*r;
sensitivityAlpha(wi) = l'*data.B*r;
cL6(wi) = alpha*sensitivityAlpha(wi);
if abs(lambda(wi))>eps
    radialSensitivityAlpha(wi) = real(conj(lambda(wi))/abs(lambda(wi))*sensitivityAlpha(wi));
end
selectedRight{wi} = r;
selectedLeft{wi} = l;
matchScore(wi) = score;
fprintf('Section 1 branch %d/%d: w %.4f lambda %.8f%+.8fi.\n', ...
    wi,numel(wGrid),w,real(lambda(wi)),imag(lambda(wi)));
end

function local_plot_branch(t, directBoundary, pencilThreshold, outputDir)
fig = figure('Visible','off','Color','w','Position',[80 80 1250 900]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(t.w, t.lambdaReal, 'k-', 'LineWidth', 1.7); hold on;
yline(1, 'r--', 'boundary');
xline(directBoundary.W, 'b--', 'direct');
xline(pencilThreshold.W, 'm:', 'pencil');
xlabel('w'); ylabel('Re(\lambda_+)'); grid on; title('Tracked critical branch');
nexttile;
plot(t.w, t.cBaseReal, 'LineWidth', 1.5); hold on;
plot(t.w, t.cL6Real, 'LineWidth', 1.5);
plot(t.w, t.lambdaReal, 'k--', 'LineWidth', 1.2);
xlabel('w'); ylabel('real contribution');
legend({'Re(l^*Ar)','Re((1-w)l^*Br)','Re(\lambda)'}, 'Location','best');
grid on; title('Biorthogonal decomposition');
nexttile;
plot(t.w, t.dLambdaDAlphaReal, 'LineWidth', 1.5); hold on;
plot(t.w, -t.finiteDifferenceDwReal, '--', 'LineWidth', 1.5);
xlabel('w'); ylabel('d Re(\lambda)/d\alpha');
legend({'l^*Br','finite difference'}, 'Location','best'); grid on;
title('Causal L6 sensitivity');
nexttile;
semilogy(t.w, t.conditionNumber, 'LineWidth', 1.4); hold on;
semilogy(t.w, t.rightResidual, '--', 'LineWidth', 1.2);
semilogy(t.w, t.leftResidual, '--', 'LineWidth', 1.2);
xlabel('w'); ylabel('value'); grid on;
legend({'condition number','right residual','left residual'}, 'Location','best');
title('Numerical conditioning');
sgtitle('Section 1: continuously tracked positive branch');
l6ns_save_figure(fig, outputDir, 'critical_branch_biorthogonal');
close(fig);
end

function local_plot_asymptotic(t, edgeB, outputDir)
fig = figure('Visible','off','Color','w','Position',[100 100 1100 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(1./t.alpha, t.leadingRealOverAlpha, '-o', 'LineWidth', 1.5); hold on;
yline(edgeB, 'k--', 'positive edge of B');
xlabel('1/(1-w)'); ylabel('Re(\lambda_+)/(1-w)'); grid on;
title('Negative-w asymptotic scaling');
nexttile;
plot(t.w, t.medianPrincipalCosine, '-o', 'LineWidth', 1.5); hold on;
plot(t.w, t.minimumPrincipalCosine, '-o', 'LineWidth', 1.2);
ylim([0 1.02]); xlabel('w'); ylabel('principal-angle cosine');
legend({'median','minimum'}, 'Location','best'); grid on;
title('Positive subspace inheritance from B');
l6ns_save_figure(fig, outputDir, 'negative_w_positive_edge_scaling');
close(fig);
end

function local_plot_actions(r, actionA, actionB, uA, uB, data, outputDir)
n = data.PopulationSize;
mapSize = data.MapSize;
vectors = {r, actionA, actionB, uA, uB};
rowLabels = {'r','Ar','(1-w)Br','u_A','u_B'};
populationLabels = {'S','C','I','E'};
fig = figure('Visible','off','Color','w','Position',[80 80 1500 1350]);
tiledlayout(numel(vectors),4,'TileSpacing','compact','Padding','compact');
for row = 1:numel(vectors)
    v = real(vectors{row});
    maps = {reshape(v(1:n),mapSize), reshape(v(n+(1:n)),mapSize), ...
        reshape(v(2*n+(1:n)),mapSize)};
    maps{4} = (1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit = max(cellfun(@(x) max(abs(x(:))), maps));
    for col = 1:4
        nexttile;
        imagesc(maps{col}); axis image off;
        if limit > 0; clim([-limit limit]); end
        colormap(gca, parula); colorbar;
        if row == 1; title(populationLabels{col}); end
        if col == 1; ylabel(rowLabels{row}, 'FontWeight','bold'); end
    end
end
sgtitle('Critical mode: base and L6 vector actions');
l6ns_save_figure(fig, outputDir, 'critical_mode_vector_actions');
close(fig);
end

function local_plot_ablation(t, outputDir)
fig = figure('Visible','off','Color','w','Position',[100 100 900 500]);
bar(categorical(t.removedBlock), t.wCriticalShift);
yline(0, 'k--'); ylabel('\Delta w_c after derivative-block removal');
title('Compensated L6 derivative ablations'); grid on;
l6ns_save_figure(fig, outputDir, 'critical_block_ablation_thresholds');
close(fig);
end
