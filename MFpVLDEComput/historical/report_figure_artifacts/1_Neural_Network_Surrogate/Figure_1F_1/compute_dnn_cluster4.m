function compute_dnn_cluster4(sourceFile, outputFile)
% Compute the four-leading-mode E-population envelope for the h96 DNN Jacobian.
loaded = load(sourceFile, 'Section4');
J = sparse(loaded.Section4.A);
clusterModeCount = 4;
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
for modeIndex = 1:clusterModeCount
    vectors(:, modeIndex) = vectors(:, modeIndex) / ...
        max(norm(vectors(:, modeIndex)), eps);
end
[clusterBasis, ~] = qr(vectors, 0);

populationSize = size(J, 1) / 3;
mapSide = round(sqrt(populationSize));
wC = 0.3077;
wS = 1 - wC;
clusterE = wS * clusterBasis(1:populationSize, :) + ...
    wC * clusterBasis(populationSize + (1:populationSize), :);
topEigenclusterE = reshape( ...
    sqrt(sum(abs(clusterE).^2, 2)), mapSide, mapSide);
save(outputFile, 'topEigenclusterE', 'clusterModeCount', 'eigValues', '-v7');
fprintf('Saved four-mode DNN eigencluster to %s\n', outputFile);
end
