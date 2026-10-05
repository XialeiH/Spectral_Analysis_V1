function summary = run_spatial_pathway_stability_analysis(sourceFile, figureDir, resultDir)
% Spatially resolve L6 and inhibitory contributions to one critical mode.
%
% Figures are exported only as PDF. Numerical results are kept separately.

if nargin < 1 || isempty(sourceFile)
    sourceFile = getenv('SPATIAL_PATHWAY_SOURCE');
end
if nargin < 2 || isempty(figureDir)
    figureDir = getenv('SPATIAL_PATHWAY_FIGURE_DIR');
end
if nargin < 3 || isempty(resultDir)
    resultDir = getenv('SPATIAL_PATHWAY_RESULT_DIR');
end
assert(isfile(sourceFile), 'Source MAT file does not exist: %s', sourceFile);
if ~exist(figureDir, 'dir'); mkdir(figureDir); end
if ~exist(resultDir, 'dir'); mkdir(resultDir); end

loaded = load(sourceFile, 'result');
assert(isfield(loaded, 'result') && isfield(loaded.result, 'Pathway'), ...
    'Source must contain result.Pathway.');
pathway = loaded.result.Pathway;
J = sparse(pathway.JBaseline);
J6 = sparse(pathway.J6);
JI = sparse(pathway.JI);
n = pathway.PopulationSize;
side = round(sqrt(n));
assert(side * side == n, 'Population size %d is not a square spatial map.', n);
labels = {'S', 'C', 'I'};
indices = {1:n, n + (1:n), 2*n + (1:n)};

[lambda, x, y] = critical_pair(J, 10);
[x, y] = orient_pair(x, y);
assert(abs(y' * x - 1) < 1e-8, 'Critical pair is not biorthogonally normalized.');

a6 = J6 * x;
aI = JI * x;
q6 = real(conj(y) .* a6);
qI = real(conj(y) .* aI);
s6 = sum(q6);
sI = sum(qI);
assert(s6 > 0 && sI < 0, ...
    'Expected destabilizing L6 and stabilizing inhibition, got s6=%g, sI=%g.', s6, sI);
cStab = -s6 / sI;
qResidual = q6 + cStab * qI;

q6Spatial = q6(indices{1}) + q6(indices{2}) + q6(indices{3});
qISpatial = qI(indices{1}) + qI(indices{2}) + qI(indices{3});
qResidualSpatial = q6Spatial + cStab*qISpatial;

metrics = metric_rows(q6, qI, cStab, indices, labels);
totalMetrics = map_metrics(q6Spatial, qISpatial, cStab);
metrics = [metrics; struct2table(add_label(totalMetrics, 'All populations'))];

q6Maps = split_maps(q6, indices, side);
qIMaps = split_maps(qI, indices, side);
qResidualMaps = split_maps(qResidual, indices, side);
a6Maps = split_maps(real(a6), indices, side);
aIMaps = split_maps(real(aI), indices, side);
xMaps = split_maps(real(x), indices, side);
yMaps = split_maps(real(y), indices, side);

plot_population_pair(q6Maps, scale_maps(qIMaps, -cStab), labels, ...
    'L6 contribution  q_6^p(r)', ...
    sprintf('Compensated inhibition  -%.4f q_I^p(r)', cStab), ...
    sprintf(['Spatial stability contributions to the baseline critical mode  ', ...
    '(angle 0.00 deg, contrast 100; Re lambda = %.6f)'], real(lambda)), ...
    fullfile(figureDir, '01_Population_Resolved_L6_and_Compensated_Inhibition_Stability_Maps.pdf'));

plot_total_maps(reshape(q6Spatial, side, side), ...
    reshape(-cStab*qISpatial, side, side), ...
    reshape(qResidualSpatial, side, side), cStab, real(lambda), ...
    fullfile(figureDir, '02_Population_Summed_Stability_Contribution_and_Residual_Maps.pdf'));

plot_residual_maps(qResidualMaps, labels, metrics, cStab, ...
    fullfile(figureDir, '03_Population_Resolved_Spatial_Compensation_Residuals.pdf'));

plot_scatter(q6, qI, q6Spatial, qISpatial, cStab, indices, labels, ...
    metrics, totalMetrics, ...
    fullfile(figureDir, '04_Pixelwise_L6_vs_Compensated_Inhibition_Opposition_Scatter.pdf'));

plot_population_pair(a6Maps, aIMaps, labels, ...
    'Re(J_6 x_*) pathway action', 'Re(J_I x_*) pathway action', ...
    'How each pathway transforms the same baseline critical perturbation', ...
    fullfile(figureDir, '05_Pathway_Action_on_the_Baseline_Critical_Mode.pdf'));

plot_population_pair(xMaps, yMaps, labels, ...
    'Re(x_*) right mode: perturbation pattern', ...
    'Re(y_*) left mode: receptive stability weighting', ...
    sprintf('Baseline biorthogonal critical mode  (y_*^*x_*=1, Re lambda=%.6f)', real(lambda)), ...
    fullfile(figureDir, '06_Baseline_Critical_Right_and_Left_Mode_Maps.pdf'));

fourier = fourier_analysis(reshape(q6Spatial, side, side), ...
    reshape(qISpatial, side, side), cStab);
plot_fourier(fourier, ...
    fullfile(figureDir, '07_Spatial_Frequency_Amplitude_and_Phase_Opposition.pdf'));
plot_radial(fourier, ...
    fullfile(figureDir, '08_Radial_Spatial_Frequency_Power_Comparison.pdf'));

modeComparison = matched_destabilizing_modes(pathway, cStab, n, side);
plot_mode_comparison(modeComparison, labels, ...
    fullfile(figureDir, '09_Matched_Destabilization_L6_Increase_vs_Inhibition_Decrease_Modes.pdf'));
plot_subspace_comparison(modeComparison, ...
    fullfile(figureDir, '10_Matched_Destabilization_Critical_Subspace_Principal_Angles.pdf'));

summary = table(real(lambda), imag(lambda), s6, sI, cStab, 1.302437, ...
    abs(cStab - 1.302437), totalMetrics.cMap, totalMetrics.spatialResidual, ...
    totalMetrics.oppositionCosine, totalMetrics.centeredOppositionCorrelation, ...
    totalMetrics.oppositeSignFraction, fourier.coherence, ...
    fourier.dominantPhaseDifference, modeComparison.gamma6, modeComparison.gammaI, ...
    modeComparison.alphaL6, modeComparison.alphaI, modeComparison.rightOverlap, ...
    modeComparison.leftOverlap, modeComparison.minimumRightSubspaceCosine, ...
    modeComparison.minimumLeftSubspaceCosine, ...
    'VariableNames', {'lambdaReal','lambdaImag','s6','sI','cStabComputed', ...
    'cStabReference','cStabAbsoluteDifference','cMap','spatialResidual', ...
    'oppositionCosine','centeredOppositionCorrelation','oppositeSignFraction', ...
    'fourierCoherence','dominantFourierPhaseDifferenceRad','gamma6Increase', ...
    'gammaIDecrease','matchedAlphaL6','matchedAlphaI','rightModeOverlap', ...
    'leftModeOverlap','minimumRightSubspaceCosine','minimumLeftSubspaceCosine'});

plot_summary_dashboard(summary, metrics, ...
    fullfile(figureDir, '11_Spatial_Opposition_and_Mode_Sharing_Quantitative_Summary.pdf'));

writetable(metrics, fullfile(resultDir, 'population_spatial_opposition_metrics.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
writetable(summary, fullfile(resultDir, 'spatial_pathway_stability_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(resultDir, 'spatial_pathway_stability_result.mat'), ...
    'summary', 'metrics', 'lambda', 'x', 'y', 'a6', 'aI', 'q6', 'qI', ...
    'qResidual', 'q6Spatial', 'qISpatial', 'qResidualSpatial', ...
    'cStab', 'fourier', 'modeComparison', '-v7.3');

pdfs = dir(fullfile(figureDir, '*.pdf'));
assert(numel(pdfs) == 11, 'Expected 11 PDFs, found %d.', numel(pdfs));
fprintf('Spatial pathway stability analysis complete: %d PDFs in %s\n', ...
    numel(pdfs), figureDir);
disp(summary);
end

function rows = metric_rows(q6, qI, cStab, indices, labels)
rowStruct = repmat(add_label(map_metrics(q6(indices{1}), qI(indices{1}), cStab), ...
    labels{1}), numel(indices), 1);
for p = 1:numel(indices)
    m = map_metrics(q6(indices{p}), qI(indices{p}), cStab);
    rowStruct(p) = add_label(m, labels{p});
end
rows = struct2table(rowStruct);
end

function output = add_label(input, label)
output = struct('population', string(label), 's6', input.s6, 'sI', input.sI, ...
    'cMap', input.cMap, 'spatialResidual', input.spatialResidual, ...
    'oppositionCosine', input.oppositionCosine, ...
    'centeredOppositionCorrelation', input.centeredOppositionCorrelation, ...
    'oppositeSignFraction', input.oppositeSignFraction);
end

function m = map_metrics(q6, qI, cStab)
q6 = real(q6(:)); qI = real(qI(:));
denominator = norm(q6) + cStab * norm(qI);
m = struct();
m.s6 = sum(q6);
m.sI = sum(qI);
m.cMap = -dot(qI, q6) / max(dot(qI, qI), eps);
m.spatialResidual = norm(q6 + cStab*qI) / max(denominator, eps);
m.oppositionCosine = dot(q6, -qI) / max(norm(q6)*norm(qI), eps);
if std(q6) > 0 && std(qI) > 0
    correlation = corrcoef(q6, -qI);
    m.centeredOppositionCorrelation = correlation(1, 2);
else
    m.centeredOppositionCorrelation = NaN;
end
weights = abs(q6) + cStab*abs(qI);
m.oppositeSignFraction = sum(weights .* (q6.*qI < 0)) / max(sum(weights), eps);
end

function maps = split_maps(vector, indices, side)
maps = cell(1, numel(indices));
for p = 1:numel(indices)
    maps{p} = reshape(vector(indices{p}), side, side);
end
end

function output = scale_maps(input, scalar)
output = cell(size(input));
for k = 1:numel(input); output{k} = scalar * input{k}; end
end

function [lambda, x, y] = critical_pair(matrix, poolSize)
dimension = size(matrix, 1);
poolSize = min(poolSize, dimension - 2);
opts = struct('tol', 1e-10, 'maxit', 2400, 'p', min(dimension, max(60, 3*poolSize)));
[right, dRight] = eigs(matrix, poolSize, 'largestreal', opts);
lambdaPool = diag(dRight);
[~, index] = max(real(lambdaPool));
lambda = lambdaPool(index);
x = right(:, index);
[leftPool, dLeft] = eigs(matrix', poolSize, 'largestreal', opts);
leftLambda = diag(dLeft);
[~, leftIndex] = min(abs(leftLambda - conj(lambda)));
y = leftPool(:, leftIndex);
x = x / norm(x);
beta = y' * x;
assert(abs(beta) > 1e-10, 'Critical left/right eigenvectors are nearly orthogonal.');
y = y / conj(beta);
end

function [x, y] = orient_pair(x, y)
[~, anchor] = max(abs(x));
phase = angle(x(anchor));
x = x * exp(-1i * phase);
y = y * exp(-1i * phase);
end

function plot_population_pair(topMaps, bottomMaps, labels, topLabel, bottomLabel, overallTitle, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [40 40 1500 850]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:3
    limit = max(abs([topMaps{p}(:); bottomMaps{p}(:)]));
    limit = max(limit, eps);
    ax = nexttile(t, p); imagesc(ax, topMaps{p}); axis(ax, 'image');
    colormap(ax, parula(256)); clim(ax, [-limit limit]); colorbar(ax);
    title(ax, sprintf('%s target: %s', labels{p}, topLabel), 'FontWeight', 'bold');
    xlabel(ax, 'spatial column'); ylabel(ax, 'spatial row');
    ax = nexttile(t, 3+p); imagesc(ax, bottomMaps{p}); axis(ax, 'image');
    colormap(ax, parula(256)); clim(ax, [-limit limit]); colorbar(ax);
    title(ax, sprintf('%s target: %s', labels{p}, bottomLabel), 'FontWeight', 'bold');
    xlabel(ax, 'spatial column'); ylabel(ax, 'spatial row');
end
title(t, overallTitle, 'FontSize', 15, 'FontWeight', 'bold');
export_pdf(fig, file); close(fig);
end

function plot_total_maps(q6, compensatedI, residual, cStab, lambdaReal, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 470]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
maps = {q6, compensatedI, residual};
names = {'q_6(r): L6 destabilization', ...
    sprintf('-%.4f q_I(r): compensated inhibition', cStab), ...
    'q_{res}(r)=q_6(r)+c_{stab}q_I(r)'};
limit = max(abs([q6(:); compensatedI(:)])); limit = max(limit, eps);
for k = 1:3
    ax = nexttile(t); imagesc(ax, maps{k}); axis(ax, 'image');
    colormap(ax, parula(256)); clim(ax, [-limit limit]); colorbar(ax);
    title(ax, names{k}, 'FontWeight', 'bold');
    xlabel(ax, 'spatial column'); ylabel(ax, 'spatial row');
end
title(t, sprintf('Population-summed stability maps, baseline Re lambda = %.6f', lambdaReal), ...
    'FontSize', 15, 'FontWeight', 'bold');
export_pdf(fig, file); close(fig);
end

function plot_residual_maps(maps, labels, metrics, cStab, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 470]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:3
    ax = nexttile(t); imagesc(ax, maps{p}); axis(ax, 'image');
    limit = max(abs(maps{p}(:))); limit = max(limit, eps);
    colormap(ax, parula(256)); clim(ax, [-limit limit]); colorbar(ax);
    title(ax, sprintf('%s residual, R_{spatial}^{%s}=%.3f', ...
        labels{p}, labels{p}, metrics.spatialResidual(p)), 'FontWeight', 'bold');
    xlabel(ax, 'spatial column'); ylabel(ax, 'spatial row');
end
title(t, sprintf('Population-resolved residual q_6^p + %.4f q_I^p', cStab), ...
    'FontSize', 15, 'FontWeight', 'bold');
export_pdf(fig, file); close(fig);
end

function plot_scatter(q6, qI, q6Spatial, qISpatial, cStab, indices, labels, ...
    metrics, totalMetrics, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1450 1050]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(3);
for p = 1:3
    ax = nexttile(t); x = q6(indices{p}); y = -cStab*qI(indices{p});
    scatter(ax, x, y, 13, colors(p,:), 'filled', 'MarkerFaceAlpha', 0.55); hold(ax, 'on');
    range = max(abs([x(:); y(:)])); plot(ax, [-range range], [-range range], 'k--');
    grid(ax, 'on'); axis(ax, 'square'); xlim(ax, [-range range]); ylim(ax, [-range range]);
    xlabel(ax, 'q_6^p(r)'); ylabel(ax, '-c_{stab}q_I^p(r)');
    title(ax, sprintf('%s: cosine %.3f, opposite-sign %.1f%%', labels{p}, ...
        metrics.oppositionCosine(p), 100*metrics.oppositeSignFraction(p)), 'FontWeight', 'bold');
end
ax = nexttile(t); hold(ax, 'on');
scatter(ax, q6Spatial, -cStab*qISpatial, 13, [0.2 0.2 0.2], ...
    'filled', 'MarkerFaceAlpha', 0.5);
range = max(abs([q6Spatial(:); cStab*qISpatial(:)]));
plot(ax, [-range range], [-range range], 'k--');
grid(ax, 'on'); axis(ax, 'square'); xlim(ax, [-range range]); ylim(ax, [-range range]);
ax.XAxis.Exponent = 0; ax.YAxis.Exponent = 0;
xlabel(ax, 'population-summed q_6(r)'); ylabel(ax, 'population-summed -c_{stab}q_I(r)');
title(ax, sprintf('Population-summed map: cosine %.3f, centered corr %.3f', ...
    totalMetrics.oppositionCosine, totalMetrics.centeredOppositionCorrelation), 'FontWeight', 'bold');
title(t, 'Pixelwise spatial opposition: equality lies on the dashed diagonal', ...
    'FontSize', 15, 'FontWeight', 'bold');
export_pdf(fig, file); close(fig);
end

function f = fourier_analysis(q6, qI, cStab)
f6 = fftshift(fft2(q6));
fI = fftshift(fft2(-cStab*qI));
amplitude6 = abs(f6); amplitudeI = abs(fI);
phaseDifference = angle(f6 .* conj(fI));
mask = (amplitude6 + amplitudeI) > 0.05*max(amplitude6(:) + amplitudeI(:));
weights = amplitude6 .* amplitudeI .* mask;
f = struct('F6', f6, 'FI', fI, 'Amplitude6', amplitude6, ...
    'AmplitudeI', amplitudeI, 'PhaseDifference', phaseDifference);
f.coherence = abs(sum(f6(:).*conj(fI(:)))) / ...
    max(norm(f6(:))*norm(fI(:)), eps);
[~, dominant] = max(weights(:));
f.dominantPhaseDifference = phaseDifference(dominant);
[f.RadialFrequency, f.RadialPower6] = radial_power(amplitude6.^2);
[~, f.RadialPowerI] = radial_power(amplitudeI.^2);
end

function [radius, power] = radial_power(powerMap)
[rows, columns] = size(powerMap);
[x, y] = meshgrid((1:columns)-ceil((columns+1)/2), (1:rows)-ceil((rows+1)/2));
r = round(sqrt(x.^2 + y.^2));
maximum = floor(min(rows, columns)/2);
radius = (0:maximum)'; power = nan(size(radius));
for k = 0:maximum
    values = powerMap(r == k);
    power(k+1) = mean(values);
end
power = power / max(sum(power), eps);
end

function plot_fourier(f, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 470]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(t); imagesc(ax, log10(f.Amplitude6 + eps)); axis(ax, 'image'); colorbar(ax);
title(ax, 'log_{10}|FFT(q_6)|', 'FontWeight', 'bold'); xlabel(ax, 'k_x'); ylabel(ax, 'k_y');
ax = nexttile(t); imagesc(ax, log10(f.AmplitudeI + eps)); axis(ax, 'image'); colorbar(ax);
title(ax, 'log_{10}|FFT(-c_{stab}q_I)|', 'FontWeight', 'bold'); xlabel(ax, 'k_x'); ylabel(ax, 'k_y');
ax = nexttile(t); imagesc(ax, f.PhaseDifference); axis(ax, 'image'); colorbar(ax);
colormap(ax, parula(256)); clim(ax, [-pi pi]);
title(ax, 'phase(q_6)-phase(-c_{stab}q_I)', 'FontWeight', 'bold'); xlabel(ax, 'k_x'); ylabel(ax, 'k_y');
title(t, sprintf('Fourier comparison: amplitude coherence %.3f, dominant phase difference %.3f rad', ...
    f.coherence, f.dominantPhaseDifference), 'FontSize', 15, 'FontWeight', 'bold');
export_pdf(fig, file); close(fig);
end

function plot_radial(f, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 620]);
semilogy(f.RadialFrequency, f.RadialPower6, '-o', 'LineWidth', 1.8); hold on;
semilogy(f.RadialFrequency, f.RadialPowerI, '--s', 'LineWidth', 1.8);
grid on; xlabel('radial spatial-frequency bin'); ylabel('normalized radial power');
title('Radial power spectra of the opposing stability-contribution maps', 'FontWeight', 'bold');
legend({'q_6', '-c_{stab}q_I'}, 'Location', 'best');
export_pdf(fig, file); close(fig);
end

function comparison = matched_destabilizing_modes(pathway, cStab, n, side)
delta6 = 0.05;
gamma6 = 1 + delta6;
j6 = pathway.JRest + gamma6*pathway.J6 + pathway.JI;
alphaTarget = max_real(j6);
objective = @(gammaI) max_real(pathway.JRest + pathway.J6 + gammaI*pathway.JI) - alphaTarget;
linearGuess = 1 - cStab*delta6;
lower = max(0, linearGuess - 0.2); upper = min(1, linearGuess + 0.2);
if objective(lower)*objective(upper) <= 0
    gammaI = fzero(objective, [lower upper]);
else
    gammaI = linearGuess;
end
jI = pathway.JRest + pathway.J6 + gammaI*pathway.JI;
[~, xBaseline, yBaseline] = critical_pair(pathway.JBaseline, 10);
[xBaseline, yBaseline] = orient_pair(xBaseline, yBaseline);
[lambda6, x6, y6] = track_critical_branch( ...
    pathway, 'L6', gamma6, xBaseline, yBaseline);
[lambdaI, xI, yI] = track_critical_branch( ...
    pathway, 'I', gammaI, xBaseline, yBaseline);

weights = population_balance_weights(x6, n);
rightOverlap = abs((weights.*x6)'*(weights.*xI)) / ...
    max(norm(weights.*x6)*norm(weights.*xI), eps);
leftOverlap = abs((y6./weights)'*(yI./weights)) / ...
    max(norm(y6./weights)*norm(yI./weights), eps);

[right6, left6] = leading_subspaces(j6, 4);
[rightI, leftI] = leading_subspaces(jI, 4);
rightCosines = svd(orth(right6)'*orth(rightI));
leftCosines = svd(orth(left6)'*orth(leftI));

indices = {1:n, n+(1:n), 2*n+(1:n)};
comparison = struct('gamma6', gamma6, 'gammaI', gammaI, ...
    'lambdaL6', lambda6, 'lambdaI', lambdaI, 'alphaL6', max_real(j6), ...
    'alphaI', max_real(jI), 'xL6', x6, 'xI', xI, 'yL6', y6, 'yI', yI, ...
    'rightOverlap', rightOverlap, 'leftOverlap', leftOverlap, ...
    'rightSubspaceCosines', rightCosines, 'leftSubspaceCosines', leftCosines, ...
    'minimumRightSubspaceCosine', min(rightCosines), ...
    'minimumLeftSubspaceCosine', min(leftCosines), ...
    'xL6Maps', {split_maps(real(x6), indices, side)}, ...
    'xIMaps', {split_maps(real(xI), indices, side)}, ...
    'yL6Maps', {split_maps(real(y6), indices, side)}, ...
    'yIMaps', {split_maps(real(yI), indices, side)});
end

function [lambda, x, y] = track_critical_branch(pathway, axisName, gammaEnd, x, y)
gammaGrid = linspace(1, gammaEnd, 21);
lambda = NaN;
for step = 2:numel(gammaGrid)
    gamma = gammaGrid(step);
    if strcmp(axisName, 'L6')
        matrix = pathway.JRest + gamma*pathway.J6 + pathway.JI;
    else
        matrix = pathway.JRest + pathway.J6 + gamma*pathway.JI;
    end
    [lambdaPool, rightPool, leftPool] = paired_eigen_pool(matrix, 12);
    rightCosine = abs(y' * rightPool) ./ max(norm(y)*vecnorm(rightPool), eps);
    leftCosine = abs(leftPool' * x).' ./ max(vecnorm(leftPool)*norm(x), eps);
    score = sqrt(rightCosine .* leftCosine);
    [~, selected] = max(score);
    lambda = lambdaPool(selected);
    xNext = rightPool(:, selected);
    yNext = leftPool(:, selected);
    phase = angle(x' * xNext);
    x = xNext * exp(-1i*phase);
    y = yNext * exp(-1i*phase);
end
end

function [lambda, right, left] = paired_eigen_pool(matrix, count)
dimension = size(matrix, 1);
opts = struct('tol', 1e-9, 'maxit', 2000, ...
    'p', min(dimension, max(70, 3*count)));
[right, dRight] = eigs(matrix, count, 'largestreal', opts);
lambda = diag(dRight);
[leftCandidates, dLeft] = eigs(matrix', count + 6, 'largestreal', opts);
leftLambda = diag(dLeft);
left = zeros(size(right));
for k = 1:count
    right(:,k) = right(:,k) / norm(right(:,k));
    [~, index] = min(abs(leftLambda - conj(lambda(k))));
    candidate = leftCandidates(:,index);
    beta = candidate' * right(:,k);
    left(:,k) = candidate / conj(beta);
end
end

function alpha = max_real(matrix)
opts = struct('tol', 1e-9, 'maxit', 1800, 'p', 50);
values = eigs(matrix, 4, 'largestreal', opts);
alpha = max(real(values));
end

function weights = population_balance_weights(vector, n)
scale = [norm(vector(1:n))/sqrt(n), norm(vector(n+(1:n)))/sqrt(n), ...
    norm(vector(2*n+(1:n)))/sqrt(n)];
scale = max(scale, 1e-12);
weights = [ones(n,1)/scale(1); ones(n,1)/scale(2); ones(n,1)/scale(3)];
end

function [right, left] = leading_subspaces(matrix, count)
opts = struct('tol', 1e-9, 'maxit', 1800, 'p', 60);
[right, d] = eigs(matrix, count, 'largestreal', opts);
lambda = diag(d);
[leftPool, dLeft] = eigs(matrix', count + 4, 'largestreal', opts);
lambdaLeft = diag(dLeft);
left = zeros(size(matrix,1), count);
for k = 1:count
    [~, index] = min(abs(lambdaLeft - conj(lambda(k))));
    left(:,k) = leftPool(:,index);
end
end

function plot_mode_comparison(c, labels, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [20 20 1500 1450]);
t = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
groups = {c.xL6Maps, c.xIMaps, c.yL6Maps, c.yIMaps};
prefix = {sprintf('L6 increased (gamma_6=%.3f), right x', c.gamma6), ...
    sprintf('inhibition reduced (gamma_I=%.4f), right x', c.gammaI), ...
    'L6 increased, left y', 'inhibition reduced, left y'};
for row = 1:4
    for p = 1:3
        ax = nexttile(t); map = groups{row}{p}; imagesc(ax, map); axis(ax, 'image');
        limit = max(abs(map(:))); limit = max(limit, eps);
        colormap(ax, parula(256)); clim(ax, [-limit limit]); colorbar(ax);
        title(ax, sprintf('%s: %s', labels{p}, prefix{row}), 'FontWeight', 'bold');
        xlabel(ax, 'spatial column'); ylabel(ax, 'spatial row');
    end
end
title(t, sprintf(['Matched destabilization: Re lambda L6+=%.6f, I-=%.6f; ', ...
    'weighted right overlap %.3f, left overlap %.3f'], ...
    c.alphaL6, c.alphaI, c.rightOverlap, c.leftOverlap), ...
    'FontSize', 15, 'FontWeight', 'bold');
export_pdf(fig, file); close(fig);
end

function plot_subspace_comparison(c, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 600]);
count = max(numel(c.rightSubspaceCosines), numel(c.leftSubspaceCosines));
values = nan(count, 2);
values(1:numel(c.rightSubspaceCosines), 1) = c.rightSubspaceCosines(:);
values(1:numel(c.leftSubspaceCosines), 2) = c.leftSubspaceCosines(:);
mode = 1:count;
bar(mode, values, 'grouped');
ylim([0 1.05]); grid on; xlabel('critical-subspace principal-angle index');
ylabel('principal-angle cosine');
title(sprintf(['Matched L6 increase vs inhibition decrease: leading critical subspaces ', ...
    '(gamma_6=%.3f, gamma_I=%.4f)'], c.gamma6, c.gammaI), 'FontWeight', 'bold');
legend({'right invariant subspace', 'left invariant subspace'}, 'Location', 'southwest');
export_pdf(fig, file); close(fig);
end

function plot_summary_dashboard(summary, metrics, file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1450 850]);
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; bar([summary.cStabComputed, summary.cMap]); grid on;
set(gca, 'XTickLabel', {'c_{stab}', 'c_{map}'}); ylabel('scaling coefficient');
title('Spectral vs best spatial cancellation scaling', 'FontWeight', 'bold');
nexttile; bar(1:height(metrics), metrics.spatialResidual); ylim([0 1]); grid on;
set(gca, 'XTick', 1:height(metrics), 'XTickLabel', cellstr(metrics.population));
ylabel('R_{spatial}'); title('Normalized spatial residual', 'FontWeight', 'bold');
nexttile; bar(1:height(metrics), ...
    [metrics.oppositionCosine, metrics.oppositeSignFraction]); ylim([-1 1]); grid on;
set(gca, 'XTick', 1:height(metrics), 'XTickLabel', cellstr(metrics.population));
ylabel('opposition metric'); legend({'anti-alignment cosine', 'weighted opposite-sign fraction'}, ...
    'Location', 'best'); title('Pixelwise spatial opposition', 'FontWeight', 'bold');
nexttile; bar([summary.rightModeOverlap, summary.leftModeOverlap, ...
    summary.minimumRightSubspaceCosine, summary.minimumLeftSubspaceCosine]); ylim([0 1]); grid on;
set(gca, 'XTickLabel', {'right mode', 'left mode', 'right subspace min', 'left subspace min'});
ylabel('overlap / cosine'); title('Matched destabilizing-mode similarity', 'FontWeight', 'bold');
title(t, sprintf(['Spatial mechanism summary: s_6=%.4g, s_I=%.4g, ', ...
    'c_{stab}=%.6f, Re lambda=%.6f'], summary.s6, summary.sI, ...
    summary.cStabComputed, summary.lambdaReal), 'FontSize', 15, 'FontWeight', 'bold');
export_pdf(fig, file); close(fig);
end

function export_pdf(fig, file)
set(fig, 'Renderer', 'painters');
exportgraphics(fig, file, 'ContentType', 'vector');
end
