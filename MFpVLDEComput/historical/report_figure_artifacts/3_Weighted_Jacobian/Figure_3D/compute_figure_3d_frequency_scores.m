function compute_figure_3d_frequency_scores(inputFile, outputFile)
% Compute radial spatial-frequency scores for every right singular vector.

loaded = load(inputFile, 'J', 'condition');
J = full(loaded.J);

[~, singularMatrix, rightSingularVectors] = svd(J, 'econ');
singularValues = diag(singularMatrix);
clear J singularMatrix

mapSize = 40;
populationSize = mapSize^2;
frequency = (-mapSize/2:(mapSize/2 - 1)) / mapSize;
[frequencyX, frequencyY] = meshgrid(frequency, frequency);
radialFrequencySquared = frequencyX.^2 + frequencyY.^2;

numModes = numel(singularValues);
kappa = zeros(numModes, 1);
for modeIndex = 1:numModes
    numerator = 0;
    denominator = 0;
    for populationIndex = 1:3
        indices = (populationIndex - 1) * populationSize + (1:populationSize);
        modeMap = reshape(rightSingularVectors(indices, modeIndex), mapSize, mapSize);
        modeSpectrum = fftshift(fft2(modeMap));
        modePower = abs(modeSpectrum).^2;
        numerator = numerator + sum(radialFrequencySquared .* modePower, 'all');
        denominator = denominator + sum(modePower, 'all');
    end
    kappa(modeIndex) = numerator / denominator;
end

condition = loaded.condition;
save(outputFile, 'singularValues', 'kappa', 'condition', 'frequency', '-v7.3');
end
