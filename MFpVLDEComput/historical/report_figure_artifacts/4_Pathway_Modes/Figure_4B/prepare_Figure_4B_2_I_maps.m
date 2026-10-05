function prepare_Figure_4B_2_I_maps
% Compute normalized I-population singular-mode maps for J_full and J_I.

artifactRoot = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(fileparts(artifactRoot)));
geometryFile = fullfile(projectRoot, 'geometry_analysis', 'results', ...
    'full_h96_geometry', ...
    'geometry_sections4_5_h96baseline_contr100_angle_0.00.mat');
outputFile = fullfile(artifactRoot, 'Figure_4B_2_singular_I_maps.mat');

loaded = load(geometryFile, 'J');
Jfull = loaded.J;
n = size(Jfull, 1) / 3;
assert(n == 1600, 'Expected three 1600-state populations.');
iRows = 2 * n + (1:n);

options.tol = 1e-10;
options.maxit = 2000;
options.disp = 0;

[leftFull, singularFull, rightFull] = svds(Jfull, 1, 'largest', options);
[leftI, singularI, rightIBlock] = svds(Jfull(:, iRows), 1, 'largest', options);

[leftFull, rightFull] = align_pair(leftFull, rightFull);
phaseI = leftI(iRows)' * rightIBlock;
if abs(phaseI) > 0
    leftI = leftI * exp(1i * angle(phaseI));
end

singularMaps.J_full.RightSingular = normalized_map(real(rightFull(iRows)));
singularMaps.J_full.LeftSingular = normalized_map(real(leftFull(iRows)));
singularMaps.J_I.RightSingular = normalized_map(real(rightIBlock));
singularMaps.J_I.LeftSingular = normalized_map(real(leftI(iRows)));

for rowName = {'J_full', 'J_I'}
    name = rowName{1};
    if dot(singularMaps.(name).RightSingular(:), ...
            singularMaps.(name).LeftSingular(:)) < 0
        singularMaps.(name).LeftSingular = ...
            -singularMaps.(name).LeftSingular;
    end
    pairMean = mean(singularMaps.(name).RightSingular(:)) + ...
        mean(singularMaps.(name).LeftSingular(:));
    if pairMean < 0
        singularMaps.(name).RightSingular = ...
            -singularMaps.(name).RightSingular;
        singularMaps.(name).LeftSingular = ...
            -singularMaps.(name).LeftSingular;
    end
end

metadata.Population = 'I';
metadata.AngleDeg = 0;
metadata.Contrast = 100;
metadata.Normalization = 'Each displayed map has unit L2 norm.';
metadata.GeometryFile = geometryFile;
metadata.JFullTopSingularValue = singularFull;
metadata.JITopSingularValue = singularI;
save(outputFile, 'singularMaps', 'metadata');
fprintf('Saved %s.\n', outputFile);
end

function [leftVector, rightVector] = align_pair(leftVector, rightVector)
phase = leftVector' * rightVector;
if abs(phase) > 0
    leftVector = leftVector * exp(1i * angle(phase));
end
end

function map = normalized_map(vector)
vector = vector(:) / max(norm(vector), eps);
map = reshape(vector, 40, 40);
end
