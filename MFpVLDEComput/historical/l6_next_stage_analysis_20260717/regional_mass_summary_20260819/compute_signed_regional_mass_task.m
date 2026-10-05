function compute_signed_regional_mass_task(taskIndex, sourceRoot, outputRoot)
% Compute signed regional masses for J_I early and J_6 late input/output maps.

if nargin < 1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin < 2 || isempty(sourceRoot)
    sourceRoot = getenv('REGIONAL_MASS_SOURCE_ROOT');
end
if nargin < 3 || isempty(outputRoot)
    outputRoot = getenv('REGIONAL_MASS_OUTPUT_ROOT');
end
if ~exist(outputRoot, 'dir'); mkdir(outputRoot); end

angles = [0 7.5 15 22.5];
contrasts = [19 42 66 100];
[contrastIndex, angleIndex] = ind2sub([4 4], taskIndex);
angleDeg = angles(angleIndex);
contrast = contrasts(contrastIndex);

dynamicFile = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW0p00_contr%d_angle_%.2f.mat', ...
    contrast, angleDeg));
frozenFile = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
    'geometry_sections4_5_h96baseline_L6eqW1p00_contr%d_angle_%.2f.mat', ...
    contrast, angleDeg));
dynamicData = load(dynamicFile, 'Section4');
frozenData = load(frozenFile, 'Section4');
jFull = sparse(dynamicData.Section4.A);
jFrozen = sparse(frozenData.Section4.A);
clear dynamicData frozenData

dimension = size(jFull, 1);
n = dimension / 3;
if dimension ~= 4800 || n ~= round(n)
    error('RegionalMass:Dimension', 'Expected a 4800-dimensional Jacobian.');
end
eRows = 1:2*n;
iRows = 2*n + (1:n);
j6 = sparse(jFull - jFrozen);
jI = sparse(dimension, dimension);
jI(:, iRows) = jFull(:, iRows);
clear jFull jFrozen

[iOutput, ~, iInput, singularFlag] = svds(jI, 1, 'largest', ...
    struct('tol', 1e-9, 'maxit', 3500, 'p', 80, 'disp', 0));
if singularFlag ~= 0
    error('RegionalMass:Svds', 'svds returned flag %d for J_I.', singularFlag);
end

j6Active = j6(eRows, eRows);
eigsOptions = struct('tol', 1e-10, 'maxit', 5000, ...
    'p', 80, 'disp', 0, 'isreal', true);
[rightCandidates, rightDiagonal, rightFlag] = ...
    eigs(j6Active, 16, 'largestreal', eigsOptions);
[leftCandidates, leftDiagonal, leftFlag] = ...
    eigs(j6Active', 16, 'largestreal', eigsOptions);
if rightFlag ~= 0 || leftFlag ~= 0
    error('RegionalMass:Eigs', ...
        'eigs flags for J_6 were right=%d, left=%d.', rightFlag, leftFlag);
end
rightValues = diag(rightDiagonal);
[~, rightOrder] = sortrows([real(rightValues), imag(rightValues)], [-1 -2]);
lambda = rightValues(rightOrder(1));
rightActive = rightCandidates(:, rightOrder(1));
leftValues = diag(leftDiagonal);
[~, leftIndex] = min(abs(leftValues - conj(lambda)));
leftActive = leftCandidates(:, leftIndex);

j6Right = [rightActive; j6(iRows, eRows) * rightActive / lambda];
j6Left = [leftActive; zeros(n, 1)];
overlap = j6Left' * j6Right;
if abs(overlap) <= 1e-12
    error('RegionalMass:Overlap', 'Leading J_6 left/right overlap is singular.');
end
j6Left = j6Left * exp(1i * angle(overlap));

masks = local_domain_masks(angleDeg, 40);
iOutputE = 0.6923 * iOutput(1:n) + 0.3077 * iOutput(n + (1:n));
j6InputE = 0.6923 * j6Left(1:n) + 0.3077 * j6Left(n + (1:n));
j6OutputE = 0.6923 * j6Right(1:n) + 0.3077 * j6Right(n + (1:n));
j6PreferredScore = sum(j6OutputE(masks.Preferred(:)));
if abs(j6PreferredScore) > 0
    phase = exp(-1i * angle(j6PreferredScore));
    j6InputE = j6InputE * phase;
    j6OutputE = j6OutputE * phase;
end
iInputMap = reshape(real(iInput(iRows)), 40, 40);
iOutputMap = reshape(real(iOutputE), 40, 40);
j6InputMap = reshape(real(j6InputE), 40, 40);
j6OutputMap = reshape(real(j6OutputE), 40, 40);

[iInputMap, iOutputMap] = local_align_pair(iInputMap, iOutputMap, masks.Preferred);
[j6InputMap, j6OutputMap] = ...
    local_align_pair(j6InputMap, j6OutputMap, masks.Preferred);

regionNames = {'Preferred', 'Surrounding', 'Orthogonal'};
maps = {iInputMap, iOutputMap, j6InputMap, j6OutputMap};
seriesNames = {'I_early_input', 'I_early_output', ...
    'L6_late_input', 'L6_late_output'};
mass = zeros(4, 3);
for seriesIndex = 1:4
    denominator = sum(abs(maps{seriesIndex}(:)));
    for regionIndex = 1:3
        mask = masks.(regionNames{regionIndex});
        mass(seriesIndex, regionIndex) = ...
            mean(maps{seriesIndex}(mask)) / max(denominator, eps);
    end
end

rows = cell(12, 6);
rowIndex = 0;
for seriesIndex = 1:4
    for regionIndex = 1:3
        rowIndex = rowIndex + 1;
        rows(rowIndex, :) = {angleDeg, contrast, string(seriesNames{seriesIndex}), ...
            string(regionNames{regionIndex}), mass(seriesIndex, regionIndex), ...
            nnz(masks.(regionNames{regionIndex}))};
    end
end
result = cell2table(rows, 'VariableNames', ...
    {'angleDeg', 'contrast', 'series', 'region', 'signedRegionalMass', 'maskPixels'});
tag = sprintf('angle%s_contrast%d', ...
    strrep(sprintf('%.2f', angleDeg), '.', 'p'), contrast);
writetable(result, fullfile(outputRoot, ['signed_regional_mass_' tag '.tsv']), ...
    'FileType', 'text', 'Delimiter', '\t');
save(fullfile(outputRoot, ['signed_regional_mass_' tag '.mat']), ...
    'result', 'mass', 'masks', 'angleDeg', 'contrast', 'lambda', '-v7');
end

function masks = local_domain_masks(angleDeg, side)
switch angleDeg
    case 0
        centers = [5 10; 5 7; 5 1];
    case 7.5
        centers = [4 10; 5 7; 6 1];
    case 15
        centers = [3 10; 5 5; 10 1];
    case 22.5
        centers = [1 10; 5 5; 10 1];
    otherwise
        error('RegionalMass:Angle', 'Unsupported angle %.2f.', angleDeg);
end
names = {'Preferred', 'Surrounding', 'Orthogonal'};
masks = struct();
preferred = false(side, side);
preferred(centers(1, 1), centers(1, 2)) = true;
surrounding = false(side, side);
surroundingRows = centers(2, 1) + [0 -1 1 0 0];
surroundingColumns = centers(2, 2) + [0 0 0 -1 1];
surrounding(sub2ind([side side], surroundingRows, surroundingColumns)) = true;
orthogonal = false(side, side);
orthogonal(centers(3, 1), centers(3, 2)) = true;
masks.(names{1}) = preferred;
masks.(names{2}) = surrounding;
masks.(names{3}) = orthogonal;
end

function [inputMap, outputMap] = local_align_pair(inputMap, outputMap, preferredMask)
score = sum(outputMap(preferredMask));
if abs(score) <= 100 * eps(max(norm(outputMap(:)), 1))
    score = sum(inputMap(preferredMask));
end
if score < 0
    inputMap = -inputMap;
    outputMap = -outputMap;
end
end
