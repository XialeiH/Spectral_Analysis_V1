function compute_Figure_2A_1_spectrum(sourceFile, outputDir, pathwayName, weight)
% Compute one full pathway-weighted eigenspectrum for Figure 2A.1.

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

loaded = load(sourceFile, 'result');
pathway = loaded.result.Pathway;
JRest = sparse(pathway.JRest);
J6 = sparse(pathway.J6);
JI = sparse(pathway.JI);

switch pathwayName
    case 'L6'
        J = JRest + (1 - weight) * J6 + JI;
    case 'I'
        J = JRest + J6 + (1 - weight) * JI;
    otherwise
        error('Unknown pathway name: %s', pathwayName);
end

lambda = eig(full(J), 'vector');
outputFile = fullfile(outputDir, sprintf('%s_weight_%+.1f_eigenvalues.mat', pathwayName, weight));
save(outputFile, 'lambda', 'weight', 'pathwayName', '-v7.3');
fprintf('Saved %d eigenvalues to %s\n', numel(lambda), outputFile);
end
