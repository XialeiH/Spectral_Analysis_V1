function diagnostics = diagnose_ji_left_cluster_mode_counts( ...
        sourceFile, outputFile, mapOutputFile)
% Check J_I left/right cluster population norms across mode counts.

loaded = load(sourceFile, 'Section4');
jFull = sparse(loaded.Section4.A);
dimension = size(jFull, 1);
n = dimension / 3;
eRows = 1:(2*n);
iRows = 2*n + (1:n);
jI = sparse(dimension, dimension);
jI(:, iRows) = jFull(:, iRows);
active = jI(iRows, iRows);

candidateCount = 80;
[rightVectors, rightValues] = local_sorted_eigs(active, candidateCount);
[leftVectors, leftValues] = local_sorted_eigs(active', candidateCount);
currentCount = sum(real(rightValues) >= max(real(rightValues)) - 1e-3);
counts = unique(max(1, currentCount + [-4 -2 -1 0 1 2 4 8]));

rows = zeros(numel(counts), 13);
side = round(sqrt(n));
maps = struct('ModeCounts', counts, ...
    'RightE', zeros(side, side, numel(counts)), ...
    'LeftE', zeros(side, side, numel(counts)), ...
    'RightI', zeros(side, side, numel(counts)), ...
    'LeftI', zeros(side, side, numel(counts)));
for countIndex = 1:numel(counts)
    count = counts(countIndex);
    right = complex(zeros(dimension, count));
    left = complex(zeros(dimension, count));
    for modeIndex = 1:count
        right(:, modeIndex) = [ ...
            jI(eRows, iRows) * rightVectors(:, modeIndex) / ...
                rightValues(modeIndex); ...
            rightVectors(:, modeIndex)];
        left(:, modeIndex) = [zeros(2*n, 1); leftVectors(:, modeIndex)];
    end
    right = right ./ max(vecnorm(right), eps);
    left = left ./ max(vecnorm(left), eps);
    rightBasis = local_basis(right);
    leftBasis = local_basis(left);
    structuredLeftBasis = [zeros(2*n, count); ...
        local_basis(leftVectors(:, 1:count))];

    rightEnvelope = sqrt(sum(abs(rightBasis).^2, 2));
    leftEnvelope = sqrt(sum(abs(leftBasis).^2, 2));
    structuredLeftEnvelope = sqrt(sum(abs(structuredLeftBasis).^2, 2));
    principalCosines = svd(leftBasis' * structuredLeftBasis);
    projectorDistance = sqrt(max(0, 2*count - ...
        2*norm(leftBasis' * structuredLeftBasis, 'fro')^2));
    maps.RightE(:, :, countIndex) = local_e_map(rightEnvelope, n, side);
    maps.LeftE(:, :, countIndex) = local_e_map(leftEnvelope, n, side);
    maps.RightI(:, :, countIndex) = reshape( ...
        rightEnvelope(iRows), side, side);
    maps.LeftI(:, :, countIndex) = reshape( ...
        leftEnvelope(iRows), side, side);
    rows(countIndex, :) = [count, ...
        real(rightValues(count)), real(leftValues(count)), ...
        norm(rightEnvelope(eRows)), norm(rightEnvelope(iRows)), ...
        norm(leftEnvelope(eRows)), norm(leftEnvelope(iRows)), ...
        max(leftEnvelope(eRows)), ...
        norm(structuredLeftEnvelope(eRows)), ...
        min(principalCosines), projectorDistance, ...
        size(rightBasis, 2), size(leftBasis, 2)];
end

diagnostics = array2table(rows, 'VariableNames', { ...
    'modeCount', 'lastRightRealEigenvalue', 'lastLeftRealEigenvalue', ...
    'rightENorm', 'rightINorm', 'leftENorm', 'leftINorm', ...
    'leftEMaximum', 'structuredLeftENorm', ...
    'minimumPrincipalCosine', 'projectorDistance', ...
    'rightRank', 'leftRank'});
diagnostics.currentModeCount = repmat(currentCount, height(diagnostics), 1);
writetable(diagnostics, outputFile, 'FileType', 'text', 'Delimiter', '\t');
if nargin >= 3 && ~isempty(mapOutputFile)
    save(mapOutputFile, 'maps', 'diagnostics', '-v7');
end
disp(diagnostics);
end

function [vectors, values] = local_sorted_eigs(matrix, count)
dimension = size(matrix, 1);
options = struct('tol', 1e-10, 'maxit', 5000, ...
    'p', min(dimension, max(180, 2*count + 20)), ...
    'disp', 0, 'isreal', true);
[vectors, diagonal, flag] = eigs(matrix, count, 'largestreal', options);
if flag ~= 0
    error('JIModeCount:Eigs', 'eigs returned flag %d.', flag);
end
values = diag(diagonal);
[~, order] = sortrows([real(values), imag(values)], [-1 -2]);
values = values(order);
vectors = vectors(:, order);
vectors = vectors ./ max(vecnorm(vectors), eps);
end

function basis = local_basis(vectors)
[basis, triangularFactor] = qr(vectors, 0);
rankTolerance = max(size(triangularFactor)) * ...
    eps(max(norm(triangularFactor, 2), 1));
rankValue = sum(abs(diag(triangularFactor)) > rankTolerance);
basis = basis(:, 1:rankValue);
end

function map = local_e_map(envelope, n, side)
wC = 0.3077;
e = (1-wC) * envelope(1:n) + wC * envelope(n + (1:n));
map = reshape(e, side, side);
end
