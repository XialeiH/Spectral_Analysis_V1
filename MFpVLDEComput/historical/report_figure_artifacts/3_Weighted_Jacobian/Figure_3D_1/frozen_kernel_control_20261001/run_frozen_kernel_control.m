function summary = run_frozen_kernel_control
% Joint L4/L6 pooling removal, retaining every baseline local derivative.
root = fileparts(mfilename('fullpath'));
output = fullfile(root,'results');
if ~isfolder(output); mkdir(output); end
loaded = load(fullfile(root,'frozen_derivatives.mat'));
assert(loaded.validation.RelativeFrobeniusError < 1e-10);
frozen = loaded.frozen;
baseline = loaded.baseline;
alphas = [0 0.1 0.25 0.5 0.75 1];
delta = 0.05;
summary = table;
JEndpoint = frozen_kernel_jacobian(frozen,baseline,1,1);
for index = 1:numel(alphas)
    alpha = alphas(index);
    [J,~,rowErrors] = frozen_kernel_jacobian(frozen,baseline,alpha,alpha);
    expected = (1-alpha)*loaded.J0+alpha*JEndpoint;
    assemblyError = norm(J-expected,'fro')/max(norm(J,'fro'),eps);
    assert(assemblyError < 1e-12);
    fprintf('alpha %.2f: computing all %d eigenvalues.\n',alpha,size(J,1));
    timer = tic;
    if alpha == 0
        eigenvalues = loaded.baselineEigenvalues;
    else
        eigenvalues = eig(full(J),'vector');
    end
    eigenSeconds = toc(timer);
    assert(numel(eigenvalues) == 4800 && all(isfinite(eigenvalues)));
    modulusCount = sum(abs(eigenvalues) <= delta);
    realPartCount = sum(abs(real(eigenvalues)) <= delta);
    boundaryDistance = min(abs(abs(eigenvalues)-delta));
    oneSummary = table(alpha,alpha,alpha,size(J,1),delta,modulusCount, ...
        realPartCount,min(abs(eigenvalues)),max(abs(eigenvalues)), ...
        max(real(eigenvalues)),boundaryDistance,max(rowErrors),assemblyError,eigenSeconds, ...
        'VariableNames',{'Alpha','L4Alpha','L6Alpha','N','Delta', ...
        'EigenvalueModulusCount','EigenvalueRealPartCount','MinimumModulus', ...
        'SpectralRadius','MaximumRealPart','DistanceToCountBoundary', ...
        'MaximumRowSumError','AssemblyLinearityError','EigenComputationSeconds'});
    save(fullfile(output,sprintf('alpha_%03d.mat',round(100*alpha))), ...
        'J','eigenvalues','oneSummary','alpha','delta','-v7.3');
    summary = [summary; oneSummary]; %#ok<AGROW>
    writetable(summary,fullfile(output,'summary.csv'));
    fprintf('alpha %.2f: disk %d, real strip %d, eig %.1f s.\n', ...
        alpha,modulusCount,realPartCount,eigenSeconds);
end
assert(height(summary) == 6 && all(summary.N == 4800));
save(fullfile(output,'summary.mat'),'summary','alphas','delta');
disp(summary);
fprintf('Completed frozen-derivative kernel intervention; no equilibria were moved.\n');
end
