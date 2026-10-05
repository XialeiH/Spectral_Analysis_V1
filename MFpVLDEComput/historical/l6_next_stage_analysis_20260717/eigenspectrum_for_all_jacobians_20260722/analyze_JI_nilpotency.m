function analyze_JI_nilpotency(sourceRoot, outputFile)
% Audit whether the inhibition derivative component is zero or nilpotent.

angles = [0 7.5 15 22.5];
contrast = 100;
populationNames = {'S', 'C', 'I'};
stats = repmat(struct(), numel(angles), 1);

for angleIndex = 1:numel(angles)
    angleValue = angles(angleIndex);
    file0 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
        ['geometry_sections4_5_h96baseline_L6eqWS0p00_C0p00_I0p00_' ...
        'contr%d_angle_%.2f.mat'], contrast, angleValue));
    file1 = fullfile(sourceRoot, 'geometry_mat', sprintf( ...
        ['geometry_sections4_5_h96baseline_L6eqWS0p00_C0p00_I1p00_' ...
        'contr%d_angle_%.2f.mat'], contrast, angleValue));
    loaded0 = load(file0, 'Section4');
    loaded1 = load(file1, 'Section4');
    JI = sparse(loaded0.Section4.A - loaded1.Section4.A);
    clear loaded0 loaded1

    matrixSize = size(JI, 1);
    if size(JI, 2) ~= matrixSize || mod(matrixSize, 3) ~= 0
        error('JIAudit:Size', 'Expected a square three-population Jacobian.');
    end
    populationSize = matrixSize / 3;
    blockFrobenius = zeros(3);
    blockNonzeros = zeros(3);
    for rowPopulation = 1:3
        rowIndex = (rowPopulation - 1) * populationSize + (1:populationSize);
        for columnPopulation = 1:3
            columnIndex = (columnPopulation - 1) * populationSize + ...
                (1:populationSize);
            block = JI(rowIndex, columnIndex);
            blockFrobenius(rowPopulation, columnPopulation) = norm(block, 'fro');
            blockNonzeros(rowPopulation, columnPopulation) = nnz(block);
        end
    end

    squared = JI * JI;
    frobeniusNorm = norm(JI, 'fro');
    squaredFrobeniusNorm = norm(squared, 'fro');
    topSingularValues = svds(JI, 6);
    directMaxAbsEigenvalue = NaN;
    directNonzeroEigenvalues = NaN;
    if angleIndex == 1
        directEigenvalues = eig(full(JI), 'vector');
        directMaxAbsEigenvalue = max(abs(directEigenvalues));
        directNonzeroEigenvalues = nnz(directEigenvalues ~= 0);
    end

    stats(angleIndex).Angle = angleValue;
    stats(angleIndex).Contrast = contrast;
    stats(angleIndex).MatrixSize = matrixSize;
    stats(angleIndex).Nonzeros = nnz(JI);
    stats(angleIndex).MaxAbsEntry = full(max(abs(JI(:))));
    stats(angleIndex).FrobeniusNorm = frobeniusNorm;
    stats(angleIndex).SpectralNorm = topSingularValues(1);
    stats(angleIndex).TopSingularValues = topSingularValues(:).';
    stats(angleIndex).DiagonalNonzeros = nnz(diag(JI));
    stats(angleIndex).MaxAbsDiagonal = max(abs(diag(JI)));
    stats(angleIndex).SquaredNonzeros = nnz(squared);
    stats(angleIndex).SquaredFrobeniusNorm = squaredFrobeniusNorm;
    stats(angleIndex).RelativeSquaredNorm = squaredFrobeniusNorm / ...
        max(frobeniusNorm ^ 2, eps);
    stats(angleIndex).BlockFrobenius = blockFrobenius;
    stats(angleIndex).BlockNonzeros = blockNonzeros;
    stats(angleIndex).DirectMaxAbsEigenvalue = directMaxAbsEigenvalue;
    stats(angleIndex).DirectNonzeroEigenvalues = directNonzeroEigenvalues;

    fprintf('\nangle %.2f deg, contrast %d\n', angleValue, contrast);
    fprintf('nnz=%d, max|entry|=%.17g, ||J_I||_F=%.17g, ||J_I||_2=%.17g\n', ...
        stats(angleIndex).Nonzeros, stats(angleIndex).MaxAbsEntry, ...
        frobeniusNorm, stats(angleIndex).SpectralNorm);
    fprintf('nnz(J_I^2)=%d, ||J_I^2||_F=%.17g, relative=%.17g\n', ...
        stats(angleIndex).SquaredNonzeros, squaredFrobeniusNorm, ...
        stats(angleIndex).RelativeSquaredNorm);
    fprintf('block Frobenius norms, rows/columns S C I:\n');
    disp(array2table(blockFrobenius, 'VariableNames', populationNames, ...
        'RowNames', populationNames));
    fprintf('block nonzero counts, rows/columns S C I:\n');
    disp(array2table(blockNonzeros, 'VariableNames', populationNames, ...
        'RowNames', populationNames));
    if angleIndex == 1
        fprintf('direct eig: max|lambda|=%.17g, exact nonzero count=%d\n', ...
            directMaxAbsEigenvalue, directNonzeroEigenvalues);
    end
end

save(outputFile, 'stats', 'populationNames', '-v7.3');
fprintf('\nSaved J_I audit to %s\n', outputFile);
end
