function summary = run_ji_component_structure_task(taskIndex, sourceRoot, ~, outputRoot)
% Diagnose J_I formed by retaining the inhibitory-source columns of J.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('JI_SOURCE_ROOT');
end
if nargin < 4 || isempty(outputRoot)
    outputRoot = getenv('JI_OUTPUT_ROOT');
end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
if ~isscalar(taskIndex) || taskIndex ~= round(taskIndex) || ...
        taskIndex < 1 || taskIndex > numel(angles) * numel(contrasts)
    error('JIAnalysis:TaskIndex', 'Task index must be an integer from 1 to 16.');
end
[contrastIndex, angleIndex] = ind2sub([numel(contrasts), numel(angles)], taskIndex);
angleValue = angles(angleIndex);
contrastValue = contrasts(contrastIndex);
if ~isfolder(sourceRoot)
    error('JIAnalysis:SourceRoot', 'Missing source root: %s', sourceRoot);
end
if isempty(outputRoot)
    error('JIAnalysis:OutputRoot', 'JI_OUTPUT_ROOT is required.');
end
dataRoot = fullfile(outputRoot, 'data');
figureRoot = fullfile(outputRoot, 'figures');
if ~exist(dataRoot, 'dir'); mkdir(dataRoot); end
if ~exist(figureRoot, 'dir'); mkdir(figureRoot); end

jBaseline = local_baseline_jacobian(sourceRoot, angleValue, contrastValue);
dimension = size(jBaseline, 1);
if size(jBaseline, 2) ~= dimension || mod(dimension, 3) ~= 0
    error('JIAnalysis:Dimension', 'Expected a square three-population Jacobian.');
end
n = dimension / 3;
sRows = 1:n;
cRows = n + (1:n);
iRows = 2*n + (1:n);
excRows = [sRows cRows];

jI = sparse(dimension, dimension);
jI(:, iRows) = jBaseline(:, iRows);
clear jBaseline
kEI = jI(excRows, iRows);
kII = jI(iRows, iRows);
activeColumn = jI(:, iRows);
froNorm = norm(jI, 'fro');
forbiddenNonIColumnRelative = norm(jI(:, excRows), 'fro') / max(froNorm, eps);
blockRepresentation = [sparse(2*n, 2*n), kEI; sparse(n, 2*n), kII];
blockRepresentationError = norm(jI - blockRepresentation, 'fro') / max(froNorm, eps);

% The full spectrum is 2N exact zeros plus the spectrum of K_II.
[vII, lambdaII] = eig(full(kII), 'vector');
[~, realOrder] = sort(real(lambdaII), 'descend');
topRealIndex = realOrder(1);
lambdaIITopReal = lambdaII(topRealIndex);
qTopReal = local_canonical_phase(vII(:, topRealIndex));
[spectralRadiusII, absoluteIndex] = max(abs(lambdaII));
lambdaIIAbsolute = lambdaII(absoluteIndex);
qAbsolute = vII(:, absoluteIndex);
qAbsolute = local_canonical_phase(qAbsolute);
otherAbsolute = true(size(lambdaII));
otherAbsolute(absoluteIndex) = false;
absoluteNearestEigenvalueDistance = min(abs(lambdaII(otherAbsolute) - lambdaIIAbsolute));
absoluteMagnitudeGapRatio = spectralRadiusII / ...
    max(max(abs(lambdaII(otherAbsolute))), eps);
minimumAbsLambdaII = min(abs(lambdaII));
maximumRealLambdaII = max(real(lambdaII));
minimumRealLambdaII = min(real(lambdaII));
maximumImagAbsLambdaII = max(abs(imag(lambdaII)));
fullMaximumRealLambda = max(0, maximumRealLambdaII);
zeroAlgebraicMultiplicity = 2*n + sum(abs(lambdaII) <= 1e-12);
zeroGeometricMultiplicityLowerBound = 2*n;
otherTopReal = true(size(lambdaII));
otherTopReal(topRealIndex) = false;
topRealNearestEigenvalueDistance = min(abs(lambdaII(otherTopReal) - lambdaIITopReal));
topRealPartGap = real(lambdaIITopReal) - max(real(lambdaII(otherTopReal)));
topRealCountWithin1e3 = sum(real(lambdaII) >= real(lambdaIITopReal) - 1e-3);

% The small top-real mode is distinct from the exact zero eigenspace only
% when K_II supplies a nonzero eigenvalue with positive real part.
if abs(lambdaIITopReal) > 1e-12
    rTopReal = [kEI * qTopReal / lambdaIITopReal; qTopReal];
    rTopReal = rTopReal / norm(rTopReal);
    topRealEigenResidual = norm(jI * rTopReal - lambdaIITopReal * rTopReal) / ...
        max(froNorm * norm(rTopReal), eps);
    topRealEnergy = local_population_energy(rTopReal, n);
    topRealSRoughness = local_periodic_roughness(rTopReal(sRows));
    topRealCRoughness = local_periodic_roughness(rTopReal(cRows));
    topRealIRoughness = local_periodic_roughness(rTopReal(iRows));
    topRealSEffectivePixels = local_effective_pixels(rTopReal(sRows));
    topRealCEffectivePixels = local_effective_pixels(rTopReal(cRows));
    topRealIEffectivePixels = local_effective_pixels(rTopReal(iRows));
    [topRealFourierTop4Fraction, topRealFourierEffectiveModes, ...
        topRealDominantKx, topRealDominantKy] = ...
        local_fourier_summary(rTopReal(iRows));
else
    rTopReal = complex(nan(dimension, 1));
    topRealEigenResidual = NaN;
    topRealEnergy = nan(1, 3);
    topRealSRoughness = NaN;
    topRealCRoughness = NaN;
    topRealIRoughness = NaN;
    topRealSEffectivePixels = NaN;
    topRealCEffectivePixels = NaN;
    topRealIEffectivePixels = NaN;
    topRealFourierTop4Fraction = NaN;
    topRealFourierEffectiveModes = NaN;
    topRealDominantKx = NaN;
    topRealDominantKy = NaN;
end

% For lambda ~= 0, q in the I block lifts to [K_EI q/lambda; q].
rAbsolute = [kEI * qAbsolute / lambdaIIAbsolute; qAbsolute];
rAbsolute = rAbsolute / norm(rAbsolute);
absoluteEigenResidual = norm(jI * rAbsolute - lambdaIIAbsolute * rAbsolute) / ...
    max(froNorm * norm(rAbsolute), eps);
absoluteEnergy = local_population_energy(rAbsolute, n);
absoluteSRoughness = local_periodic_roughness(rAbsolute(sRows));
absoluteCRoughness = local_periodic_roughness(rAbsolute(cRows));
absoluteIRoughness = local_periodic_roughness(rAbsolute(iRows));
absoluteSEffectivePixels = local_effective_pixels(rAbsolute(sRows));
absoluteCEffectivePixels = local_effective_pixels(rAbsolute(cRows));
absoluteIEffectivePixels = local_effective_pixels(rAbsolute(iRows));
[absoluteFourierTop4Fraction, absoluteFourierEffectiveModes, ...
    absoluteDominantKx, absoluteDominantKy] = ...
    local_fourier_summary(rAbsolute(iRows));

% Singular vectors measure the strongest one-step inhibitory-source action.
opts = struct('tol', 1e-10, 'maxit', 3000, 'disp', 0);
[uSingular, singularMatrix, vSingular] = svds(activeColumn, 4, 'largest', opts);
singularValues = diag(singularMatrix);
[singularValues, singularOrder] = sort(real(singularValues), 'descend');
uSingular = uSingular(:, singularOrder);
vSingular = vSingular(:, singularOrder);
sigma1 = singularValues(1);
sigma2 = singularValues(2);
sigma3 = singularValues(3);
sigma4 = singularValues(4);
singularGapRatio = sigma1 / max(sigma2, eps);
stableRank = froNorm^2 / max(sigma1^2, eps);
u1 = uSingular(:, 1);
vI1 = vSingular(:, 1);
[u1, vI1] = local_canonical_singular_sign(u1, vI1);
singularResidual = norm(activeColumn * vI1 - sigma1 * u1) / ...
    max(sigma1 * norm(u1), eps);
adjointSingularResidual = norm(activeColumn' * u1 - sigma1 * vI1) / ...
    max(sigma1 * norm(vI1), eps);
singularLeftEnergy = local_population_energy(u1, n);
singularInputRoughness = local_periodic_roughness(vI1);
singularOutputSRoughness = local_periodic_roughness(u1(sRows));
singularOutputCRoughness = local_periodic_roughness(u1(cRows));
singularOutputIRoughness = local_periodic_roughness(u1(iRows));

% Two unrelated excitatory vectors are both exact zero eigenvectors.
rng(2906 + taskIndex, 'twister');
nullVectorA = [randn(2*n, 1); zeros(n, 1)];
nullVectorB = [randn(2*n, 1); zeros(n, 1)];
nullVectorA = nullVectorA / norm(nullVectorA);
nullVectorB = nullVectorB / norm(nullVectorB);
nullResidualA = norm(jI * nullVectorA) / max(froNorm, eps);
nullResidualB = norm(jI * nullVectorB) / max(froNorm, eps);
nullVectorOverlap = abs(nullVectorA' * nullVectorB);

% Block norms expose which target populations receive the I-source pathway.
normSI = norm(jI(sRows, iRows), 'fro');
normCI = norm(jI(cRows, iRows), 'fro');
normII = norm(kII, 'fro');

summary = table(taskIndex, angleValue, contrastValue, dimension, n, nnz(jI), ...
    froNorm, norm(jI, 1), norm(jI, inf), max(abs(nonzeros(jI))), ...
    forbiddenNonIColumnRelative, blockRepresentationError, ...
    normSI, normCI, normII, ...
    zeroAlgebraicMultiplicity, zeroGeometricMultiplicityLowerBound, ...
    fullMaximumRealLambda, real(lambdaIITopReal), imag(lambdaIITopReal), ...
    topRealNearestEigenvalueDistance, topRealPartGap, topRealCountWithin1e3, ...
    topRealEigenResidual, topRealEnergy(1), topRealEnergy(2), topRealEnergy(3), ...
    topRealSRoughness, topRealCRoughness, topRealIRoughness, ...
    topRealSEffectivePixels, topRealCEffectivePixels, topRealIEffectivePixels, ...
    topRealFourierTop4Fraction, topRealFourierEffectiveModes, ...
    topRealDominantKx, topRealDominantKy, ...
    minimumRealLambdaII, maximumRealLambdaII, maximumImagAbsLambdaII, ...
    minimumAbsLambdaII, spectralRadiusII, ...
    real(lambdaIIAbsolute), imag(lambdaIIAbsolute), ...
    absoluteNearestEigenvalueDistance, absoluteMagnitudeGapRatio, ...
    absoluteEigenResidual, ...
    absoluteEnergy(1), absoluteEnergy(2), absoluteEnergy(3), ...
    absoluteSRoughness, absoluteCRoughness, absoluteIRoughness, ...
    absoluteSEffectivePixels, absoluteCEffectivePixels, absoluteIEffectivePixels, ...
    absoluteFourierTop4Fraction, absoluteFourierEffectiveModes, ...
    absoluteDominantKx, absoluteDominantKy, ...
    sigma1, sigma2, sigma3, sigma4, singularGapRatio, stableRank, ...
    singularResidual, adjointSingularResidual, ...
    singularLeftEnergy(1), singularLeftEnergy(2), singularLeftEnergy(3), ...
    singularInputRoughness, singularOutputSRoughness, ...
    singularOutputCRoughness, singularOutputIRoughness, ...
    nullResidualA, nullResidualB, nullVectorOverlap, ...
    'VariableNames', {'taskIndex','angle','contrast','dimension','populationSize', ...
    'nonzeroCount','frobeniusNorm','oneNorm','infinityNorm','maximumAbsEntry', ...
    'forbiddenNonIColumnRelative','blockRepresentationError', ...
    'frobeniusSI','frobeniusCI','frobeniusII', ...
    'zeroAlgebraicMultiplicity','zeroGeometricMultiplicityLowerBound', ...
    'fullMaximumRealLambda','topRealLambdaIIReal','topRealLambdaIIImag', ...
    'topRealNearestEigenvalueDistance','topRealPartGap','topRealCountWithin1e3', ...
    'topRealEigenResidual','topRealEnergyS','topRealEnergyC','topRealEnergyI', ...
    'topRealSRoughness','topRealCRoughness','topRealIRoughness', ...
    'topRealSEffectivePixels','topRealCEffectivePixels','topRealIEffectivePixels', ...
    'topRealFourierTop4Fraction','topRealFourierEffectiveModes', ...
    'topRealDominantKx','topRealDominantKy', ...
    'minimumRealLambdaII','maximumRealLambdaII','maximumImagAbsLambdaII', ...
    'minimumAbsLambdaII','spectralRadiusII', ...
    'absoluteLambdaIIReal','absoluteLambdaIIImag', ...
    'absoluteNearestEigenvalueDistance','absoluteMagnitudeGapRatio', ...
    'absoluteEigenResidual', ...
    'absoluteEnergyS','absoluteEnergyC','absoluteEnergyI', ...
    'absoluteSRoughness','absoluteCRoughness','absoluteIRoughness', ...
    'absoluteSEffectivePixels','absoluteCEffectivePixels','absoluteIEffectivePixels', ...
    'absoluteFourierTop4Fraction','absoluteFourierEffectiveModes', ...
    'absoluteDominantKx','absoluteDominantKy', ...
    'sigma1','sigma2','sigma3','sigma4','singularGapRatio','stableRank', ...
    'singularResidual','adjointSingularResidual', ...
    'singularLeftEnergyS','singularLeftEnergyC','singularLeftEnergyI', ...
    'singularInputRoughness','singularOutputSRoughness', ...
    'singularOutputCRoughness','singularOutputIRoughness', ...
    'nullResidualA','nullResidualB','nullVectorOverlap'});

tag = sprintf('angle%s_contrast%d', local_angle_tag(angleValue), contrastValue);
writetable(summary, fullfile(dataRoot, ['ji_structure_', tag, '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(dataRoot, ['ji_modes_', tag, '.mat']), ...
    'summary', 'lambdaII', 'lambdaIITopReal', 'rTopReal', ...
    'lambdaIIAbsolute', 'rAbsolute', ...
    'singularValues', 'u1', 'vI1', '-v7.3');
if contrastValue == 100
    local_plot_top_real_mode(rTopReal, lambdaIITopReal, ...
        topRealNearestEigenvalueDistance, topRealPartGap, angleValue, ...
        contrastValue, figureRoot, tag, n);
    local_plot_absolute_mode(rAbsolute, lambdaIIAbsolute, angleValue, ...
        contrastValue, figureRoot, tag, n);
    local_plot_singular_pair(vI1, u1, sigma1, singularGapRatio, ...
        angleValue, contrastValue, figureRoot, tag, n);
end
fprintf(['J_I angle %.2f contrast %d: zeros >=%d, maxRe=%.6g, ' ...
    'rho=%.6g, sigma1=%.6g.\n'], angleValue, contrastValue, ...
    2*n, fullMaximumRealLambda, spectralRadiusII, sigma1);
end

function local_plot_top_real_mode(vector, lambda, nearestDistance, realPartGap, ...
        angleValue, contrastValue, outputRoot, tag, n)
maps = {reshape(real(vector(1:n)), 40, 40), ...
    reshape(real(vector(n + (1:n))), 40, 40), ...
    reshape(real(vector(2*n + (1:n))), 40, 40)};
labels = {'S', 'C', 'I'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 560]);
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
colormap(fig, parula(256));
for index = 1:3
    axisHandle = nexttile(layout);
    imagesc(axisHandle, maps{index});
    axis(axisHandle, 'image');
    axis(axisHandle, 'off');
    limit = max(abs(maps{index}(:)));
    if limit <= eps; limit = 1; end
    clim(axisHandle, [-limit limit]);
    colorbar(axisHandle);
    title(axisHandle, labels{index}, 'FontSize', 12, 'FontWeight', 'bold');
end
title(layout, { ...
    sprintf('Top-real eigenmode of inhibitory-source J_I | angle %.2f deg | contrast %d', ...
        angleValue, contrastValue), ...
    sprintf('lambda=%.6g%+.6gi | nearest spectral distance %.3g | real-part gap %.3g', ...
        real(lambda), imag(lambda), nearestDistance, realPartGap)}, ...
    'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig, fullfile(outputRoot, ['ji_top_real_eigenmode_', tag, '.pdf']), ...
    'ContentType', 'vector');
close(fig);
end

function jBaseline = local_baseline_jacobian(sourceRoot, angleValue, contrastValue)
tag = 'L6eqWS0p00_C0p00_I0p00';
file = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_%s_contr%d_angle_%.2f.mat', ...
    tag, contrastValue, angleValue));
if ~isfile(file)
    error('JIAnalysis:MissingBaseline', 'Missing baseline file: %s', file);
end
loaded = load(file, 'Section4');
jBaseline = loaded.Section4.A;
end

function energy = local_population_energy(vector, n)
energy = [sum(abs(vector(1:n)).^2), ...
    sum(abs(vector(n + (1:n))).^2), ...
    sum(abs(vector(2*n + (1:n))).^2)];
energy = energy / max(sum(energy), eps);
end

function value = local_periodic_roughness(vector)
side = round(sqrt(numel(vector)));
map = reshape(real(vector), side, side);
dx = map - circshift(map, [0 1]);
dy = map - circshift(map, [1 0]);
value = (sum(dx(:).^2) + sum(dy(:).^2)) / max(sum(map(:).^2), eps);
end

function value = local_effective_pixels(vector)
energy = abs(vector(:)).^2;
energy = energy / max(sum(energy), eps);
value = 1 / max(sum(energy.^2), eps);
end

function [top4Fraction, effectiveModes, dominantKx, dominantKy] = ...
        local_fourier_summary(vector)
side = round(sqrt(numel(vector)));
map = reshape(vector, side, side);
power = abs(fftshift(fft2(map))).^2;
power = power / max(sum(power(:)), eps);
sortedPower = sort(power(:), 'descend');
top4Fraction = sum(sortedPower(1:min(4, numel(sortedPower))));
positivePower = power(power > 0);
effectiveModes = exp(-sum(positivePower .* log(positivePower)));
[~, dominantIndex] = max(power(:));
[dominantY, dominantX] = ind2sub([side side], dominantIndex);
center = side / 2 + 1;
dominantKx = dominantX - center;
dominantKy = dominantY - center;
end

function vector = local_canonical_phase(vector)
[~, index] = max(abs(vector));
phase = vector(index) / max(abs(vector(index)), eps);
vector = vector / phase;
end

function [u, v] = local_canonical_singular_sign(u, v)
[~, index] = max(abs(u));
phase = u(index) / max(abs(u(index)), eps);
u = u / phase;
v = v / phase;
end

function local_plot_absolute_mode(vector, lambda, angleValue, contrastValue, ...
        outputRoot, tag, n)
maps = {reshape(real(vector(1:n)), 40, 40), ...
    reshape(real(vector(n + (1:n))), 40, 40), ...
    reshape(real(vector(2*n + (1:n))), 40, 40)};
labels = {'S', 'C', 'I'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1450 560]);
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
colormap(fig, parula(256));
for index = 1:3
    axisHandle = nexttile(layout);
    imagesc(axisHandle, maps{index});
    axis(axisHandle, 'image');
    axis(axisHandle, 'off');
    limit = max(abs(maps{index}(:)));
    if limit <= eps; limit = 1; end
    clim(axisHandle, [-limit limit]);
    colorbar(axisHandle);
    title(axisHandle, labels{index}, 'FontSize', 12, 'FontWeight', 'bold');
end
title(layout, { ...
    sprintf('Largest-|lambda| eigenmode of inhibitory-source J_I | angle %.2f deg | contrast %d', ...
        angleValue, contrastValue), ...
    sprintf('lambda=%.6g%+.6gi | this is the strongest intrinsic K_{II} eigenmode, not the zero spectral edge', ...
        real(lambda), imag(lambda))}, ...
    'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig, fullfile(outputRoot, ['ji_largest_absolute_eigenmode_', tag, '.pdf']), ...
    'ContentType', 'vector');
close(fig);
end

function local_plot_singular_pair(vI, u, sigma1, gapRatio, angleValue, ...
        contrastValue, outputRoot, tag, n)
maps = {reshape(real(vI), 40, 40), ...
    reshape(real(u(1:n)), 40, 40), ...
    reshape(real(u(n + (1:n))), 40, 40), ...
    reshape(real(u(2*n + (1:n))), 40, 40)};
labels = {'right singular input: I', 'left singular output: S', ...
    'left singular output: C', 'left singular output: I'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1800 560]);
layout = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
colormap(fig, parula(256));
for index = 1:4
    axisHandle = nexttile(layout);
    imagesc(axisHandle, maps{index});
    axis(axisHandle, 'image');
    axis(axisHandle, 'off');
    limit = max(abs(maps{index}(:)));
    if limit <= eps; limit = 1; end
    clim(axisHandle, [-limit limit]);
    colorbar(axisHandle);
    title(axisHandle, labels{index}, 'FontSize', 12, 'FontWeight', 'bold');
end
title(layout, { ...
    sprintf('Dominant one-step singular pathway of inhibitory-source J_I | angle %.2f deg | contrast %d', ...
        angleValue, contrastValue), ...
    sprintf('sigma_1=%.6g | sigma_1/sigma_2=%.6g', sigma1, gapRatio)}, ...
    'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');
exportgraphics(fig, fullfile(outputRoot, ['ji_top_singular_pair_', tag, '.pdf']), ...
    'ContentType', 'vector');
close(fig);
end

function tag = local_angle_tag(value)
tag = strrep(sprintf('%.2f', value), '.', 'p');
end
