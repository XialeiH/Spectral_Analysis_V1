function export_angle0_sc_reduced_eigenvalues(sourceFile, outputFile)
% Export the naive 3200-state E-I eigenspectrum for angle 0, contrast 100.

data = load(sourceFile, 'Section4', 'Angle', 'Contr');
J = data.Section4.A;
n = size(J, 1) / 3;
if n ~= round(n)
    error('SCReduction:Dimension', 'Expected a 3-population Jacobian.');
end

wC = 0.3077;
wS = 1 - wC;
s = 1:n;
c = n + (1:n);
i = 2 * n + (1:n);

jEE = wS * (J(s, s) + J(s, c)) + wC * (J(c, s) + J(c, c));
jEI = wS * J(s, i) + wC * J(c, i);
jIE = J(i, s) + J(i, c);
jII = J(i, i);
JReduced = [jEE, jEI; jIE, jII];
clear J jEE jEI jIE jII

reducedEigenvalues = eig(full(JReduced), 'vector');
nearZeroThreshold = 0.05;
nearZeroCount = nnz(abs(reducedEigenvalues) <= nearZeroThreshold);
if nearZeroCount ~= 1369
    error('SCReduction:CountMismatch', ...
        'Expected 1369 reduced near-zero eigenvalues, found %d.', nearZeroCount);
end

angle = data.Angle;
contrast = data.Contr;
fullDimension = 3 * n;
reducedDimension = 2 * n;
save(outputFile, 'reducedEigenvalues', 'angle', 'contrast', ...
    'wS', 'wC', 'fullDimension', 'reducedDimension', ...
    'nearZeroThreshold', 'nearZeroCount', '-v7.3');
end
