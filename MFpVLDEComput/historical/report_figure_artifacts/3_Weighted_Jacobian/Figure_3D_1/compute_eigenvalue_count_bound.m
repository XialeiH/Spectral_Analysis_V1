function summary = compute_eigenvalue_count_bound
% Weyl prefix-product bound for the exact spectra used by Figure 3D.1.

outputDir = fileparts(mfilename('fullpath'));
artifactRoot = fileparts(outputDir);
sourceNames = ["baseline"; "remove_L4_and_L6_smoothing"];
conditionNames = ["Baseline"; "Flattened"];
frequencyNames = ["baseline"; "flattened"];
delta = 0.05;
summary = table;
curves = table;

for conditionIndex = 1:2
    sourceFile = fullfile(outputDir, 'source', sourceNames(conditionIndex)+'.mat');
    source = load(sourceFile, 'singularValues', 'eigenvalues');
    original = load(fullfile(artifactRoot, 'Figure_3D', ...
        frequencyNames(conditionIndex)+'_frequency_scores.mat'), 'singularValues');
    sigma = sort(source.singularValues(:), 'descend');
    lambdaModulus = sort(abs(source.eigenvalues(:)), 'descend');
    n = numel(sigma);
    assert(n == 4800 && numel(lambdaModulus) == n);
    assert(all(isfinite(sigma)) && all(sigma > 0));
    assert(all(isfinite(lambdaModulus)) && all(lambdaModulus > 0));
    sourceMatchError = max(abs(sigma-sort(original.singularValues(:), 'descend')));
    assert(sourceMatchError < 1e-10, 'Source spectrum must match Figure 3D.');

    k = (1:n)';
    logSigmaPrefix = cumsum(log(sigma));
    logG = logSigmaPrefix./k;
    g = exp(logG);
    logLambdaPrefix = cumsum(log(lambdaModulus));
    assert(all(logLambdaPrefix <= logSigmaPrefix+1e-8), ...
        'Computed spectra must satisfy the multiplicative Weyl inequality.');
    assert(all(lambdaModulus <= g.*(1+1e-10)));
    assert(all(diff(logG) <= 1e-12));

    firstK = find(logG <= log(delta), 1, 'first');
    if isempty(firstK)
        firstK = NaN;
        lowerBound = 0;
        gBefore = NaN;
        gAtCrossing = NaN;
    else
        lowerBound = n-firstK+1;
        gAtCrossing = g(firstK);
        gBefore = NaN;
        if firstK > 1
            gBefore = g(firstK-1);
            assert(gBefore > delta);
        end
    end
    actualCount = sum(lambdaModulus <= delta);
    singularCount = sum(sigma <= delta);
    assert(lowerBound <= actualCount);

    oneSummary = table(conditionNames(conditionIndex), n, delta, firstK, ...
        gBefore, gAtCrossing, g(end), lowerBound, actualCount, singularCount, ...
        min(lambdaModulus), sourceMatchError, string(sourceFile), ...
        'VariableNames', {'Condition','N','Delta','FirstK','GBefore', ...
        'GAtCrossing','GAtN','EigenvalueCountLowerBound','ActualEigenvalueCount', ...
        'SmallSingularValueCount','MinimumEigenvalueModulus', ...
        'SourceSingularValueMatchError','SourceFile'});
    oneCurve = table(repmat(conditionNames(conditionIndex), n, 1), k, sigma, ...
        lambdaModulus, logG, g, repmat(delta, n, 1), logG <= log(delta), ...
        'VariableNames', {'Condition','K','SingularValueDescending', ...
        'EigenvalueModulusDescending','LogG','G','Delta','GuaranteesThreshold'});
    summary = [summary; oneSummary]; %#ok<AGROW>
    curves = [curves; oneCurve]; %#ok<AGROW>
end

writetable(summary, fullfile(outputDir, 'eigenvalue_count_bound_summary.csv'));
writetable(curves, fullfile(outputDir, 'eigenvalue_count_bound_curves.csv'));
save(fullfile(outputDir, 'eigenvalue_count_bound.mat'), 'summary', 'curves', 'delta');
format long g
disp(summary(:, 1:12));
fprintf('Saved summary and all 9,600 rank-wise comparisons in %s\n', outputDir);
end
