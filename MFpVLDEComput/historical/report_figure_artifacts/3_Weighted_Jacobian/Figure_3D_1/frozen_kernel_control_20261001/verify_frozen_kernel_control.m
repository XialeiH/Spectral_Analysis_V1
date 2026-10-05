function audit = verify_frozen_kernel_control
% Independently check saved spectra, matrix identities, and local endpoint blocks.
root = fileparts(mfilename('fullpath'));
output = fullfile(root,'results');
data = load(fullfile(root,'frozen_derivatives.mat'));
source = load(data.frozen.SourceFile,'J','fixedPoint');
assert(isequal(data.frozen.FixedPoint,source.fixedPoint(:)));
assert(norm(data.J0-source.J,'fro') == 0);
summary = readtable(fullfile(output,'summary.csv'));
assert(height(summary) == 6);
audit = table;
endpoint = load(fullfile(output,'alpha_100.mat'),'J');
for index = 1:height(summary)
    alpha = summary.Alpha(index);
    one = load(fullfile(output,sprintf('alpha_%03d.mat',round(100*alpha))));
    eigenvalues = one.eigenvalues;
    assert(numel(eigenvalues) == 4800);
    assert(sum(abs(eigenvalues)<=0.05) == summary.EigenvalueModulusCount(index));
    assert(sum(abs(real(eigenvalues))<=0.05) == summary.EigenvalueRealPartCount(index));
    [rebuilt,~,rowErrors] = frozen_kernel_jacobian(data.frozen,data.baseline,alpha,alpha);
    assert(norm(rebuilt-one.J,'fro') == 0);
    expected = (1-alpha)*data.J0+alpha*endpoint.J;
    matrixError = norm(one.J-expected,'fro')/norm(one.J,'fro');
    traceValue = sum(diag(one.J));
    traceSquaredValue = sum(sum(one.J.*one.J.'));
    traceError = abs(sum(eigenvalues)-traceValue)/max(1,abs(traceValue));
    traceSquaredError = abs(sum(eigenvalues.^2)-traceSquaredValue)/max(1,abs(traceSquaredValue));
    assert(matrixError < 1e-12 && max(rowErrors) < 1e-12);
    assert(traceError < 1e-10 && traceSquaredError < 1e-10);
    audit = [audit; table(alpha,matrixError,traceError,traceSquaredError,max(rowErrors), ...
        'VariableNames',{'Alpha','MatrixLinearityError','TraceError', ...
        'TraceSquaredError','MaximumRowSumError'})]; %#ok<AGROW>
end
n = 1600;
[rows,columns] = find(endpoint.J);
assert(all(mod(rows-1,n) == mod(columns-1,n)), ...
    'Delta-kernel endpoint must have no cross-pixel coupling.');
localEigenvalues = zeros(4800,1);
for pixel = 1:n
    indices = pixel+n*(0:2);
    localEigenvalues(3*(pixel-1)+(1:3)) = eig(full(endpoint.J(indices,indices)));
end
endpointData = load(fullfile(output,'alpha_100.mat'),'eigenvalues');
endpointModulusError = max(abs(sort(abs(localEigenvalues))-sort(abs(endpointData.eigenvalues))));
assert(endpointModulusError < 1e-10);
assert(sum(abs(localEigenvalues)<=0.05) == 192);
writetable(audit,fullfile(output,'validation.csv'));
save(fullfile(output,'validation.mat'),'audit','endpointModulusError','localEigenvalues');
disp(audit);
fprintf('Verified 6 complete spectra; endpoint independent modulus error %.3g.\n',endpointModulusError);
end
