function threshold = l6ns_generalized_threshold(a, b, cfg, referenceAlpha)
% Solve (I-A)r = alpha B r through the reciprocal matrix-free pencil.

if nargin < 4
    referenceAlpha = NaN;
end
n = size(a, 1);
systemMatrix = speye(n) - a;
solver = decomposition(systemMatrix, 'lu');
operator = @(x) solver \ (b * x);

opts = struct();
opts.tol = cfg.EigsTolerance;
opts.maxit = cfg.EigsMaxIterations;
opts.p = min(max(cfg.EigsSubspaceDimension, 64), n);
opts.issym = false;
opts.isreal = isreal(a) && isreal(b);
opts.disp = 0;

[vectors, values] = eigs(operator, n, 24, 'largestreal', opts);
mu = diag(values);
alpha = 1 ./ mu;
eligible = real(alpha) > 0 & abs(imag(alpha)) < 1e-6 * max(1, abs(real(alpha)));
if ~any(eligible)
    error('No positive real generalized threshold was returned by the matrix-free pencil.');
end

indices = find(eligible);
if isfinite(referenceAlpha)
    [~, localIndex] = min(abs(real(alpha(indices)) - referenceAlpha));
else
    [~, localIndex] = min(real(alpha(indices)));
end
index = indices(localIndex);
r = vectors(:, index);
r = r / max(norm(r), eps);

threshold = struct();
threshold.Alpha = real(alpha(index));
threshold.W = 1 - threshold.Alpha;
threshold.Mu = mu(index);
threshold.RightVector = r;
threshold.Residual = norm(systemMatrix * r - threshold.Alpha * b * r) / ...
    max(norm(systemMatrix * r) + abs(threshold.Alpha) * norm(b * r), eps);
threshold.AllAlpha = alpha;
threshold.AllMu = mu;
end
