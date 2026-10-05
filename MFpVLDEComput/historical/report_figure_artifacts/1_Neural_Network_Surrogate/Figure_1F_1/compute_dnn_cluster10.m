function compute_dnn_cluster10(sourceFile, outputFile)
% Compute the ten-leading-mode E-population envelope for the h96 DNN Jacobian.
loaded = load(sourceFile, 'Section4');
J = sparse(loaded.Section4.A);
clusterModeCount = 10;
options = struct('tol', 1e-10, 'maxit', 5000, ...
    'p', 80, 'disp', 0, 'isreal', true);
[vectors, diagonal, flag] = eigs( ...
    J, clusterModeCount, 'largestreal', options);
if flag ~= 0
    error('Figure1F1:Eigs', 'eigs returned flag %d.', flag);
end
eigValues = diag(diagonal);
[~, order] = sortrows([real(eigValues), imag(eigValues)], [-1 -2]);
eigValues = eigValues(order);
vectors = vectors(:, order);
populationSize = size(J, 1) / 3;
mapSide = round(sqrt(populationSize));
wC = 0.3077;
wS = 1 - wC;
topEigenmodeE = zeros(mapSide, mapSide, clusterModeCount);
for modeIndex = 1:clusterModeCount
    vectors(:, modeIndex) = vectors(:, modeIndex) / ...
        max(norm(vectors(:, modeIndex)), eps);
    [~, pivot] = max(abs(vectors(:, modeIndex)));
    vectors(:, modeIndex) = vectors(:, modeIndex) * ...
        exp(-1i * angle(vectors(pivot, modeIndex)));
    if real(vectors(pivot, modeIndex)) < 0
        vectors(:, modeIndex) = -vectors(:, modeIndex);
    end
    eMode = wS * vectors(1:populationSize, modeIndex) + ...
        wC * vectors(populationSize + (1:populationSize), modeIndex);
    topEigenmodeE(:, :, modeIndex) = reshape(real(eMode), mapSide, mapSide);
end
[clusterBasis, ~] = qr(vectors, 0);

clusterE = wS * clusterBasis(1:populationSize, :) + ...
    wC * clusterBasis(populationSize + (1:populationSize), :);
topEigenclusterE = reshape( ...
    sqrt(sum(abs(clusterE).^2, 2)), mapSide, mapSide);
save(outputFile, 'topEigenclusterE', 'topEigenmodeE', ...
    'clusterModeCount', 'eigValues', '-v7');
fprintf('Saved ten-mode DNN eigencluster to %s\n', outputFile);
end
