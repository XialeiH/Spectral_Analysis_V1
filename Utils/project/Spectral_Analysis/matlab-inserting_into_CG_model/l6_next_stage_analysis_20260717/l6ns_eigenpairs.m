function modes = l6ns_eigenpairs(jacobian, count, selector, cfg)
% Compute matched left/right eigenpairs with biorthogonal normalization.

n = size(jacobian, 1);
count = min(count, n - 2);
opts = struct();
opts.tol = cfg.EigsTolerance;
opts.maxit = cfg.EigsMaxIterations;
opts.p = min(max(cfg.EigsSubspaceDimension, 2 * count + 8), n);
opts.disp = 0;

try
    [rightVectors, rightValues] = eigs(jacobian, count, selector, opts);
catch exception
    if n > 600
        rethrow(exception)
    end
    [rightVectors, rightValues] = local_dense_pairs(jacobian, count, selector);
end
lambda = diag(rightValues);

leftSelector = selector;
if isnumeric(selector)
    leftSelector = conj(selector);
end
try
    [leftPool, leftValues] = eigs(jacobian', count, leftSelector, opts);
catch exception
    if n > 600
        rethrow(exception)
    end
    [leftPool, leftValues] = local_dense_pairs(jacobian', count, leftSelector);
end
lambdaLeft = diag(leftValues);

used = false(numel(lambdaLeft), 1);
leftVectors = zeros(size(leftPool), 'like', leftPool);
conditionNumber = nan(count, 1);
rightResidual = nan(count, 1);
leftResidual = nan(count, 1);
jNorm = norm(jacobian, 'fro');

for k = 1:count
    distance = abs(lambdaLeft - conj(lambda(k)));
    distance(used) = inf;
    [~, index] = min(distance);
    used(index) = true;

    r = rightVectors(:, k);
    r = r / max(norm(r), eps);
    l = leftPool(:, index);
    overlap = l' * r;
    if abs(overlap) < 1e-13
        warning('Left/right overlap is nearly singular for eigenvalue %.6g%+.6gi.', ...
            real(lambda(k)), imag(lambda(k)));
    end
    l = l / conj(overlap);

    rightVectors(:, k) = r;
    leftVectors(:, k) = l;
    conditionNumber(k) = norm(l) * norm(r);
    rightResidual(k) = norm(jacobian * r - lambda(k) * r) / ...
        max((jNorm + abs(lambda(k))) * norm(r), eps);
    leftResidual(k) = norm(jacobian' * l - conj(lambda(k)) * l) / ...
        max((jNorm + abs(lambda(k))) * norm(l), eps);
end

[~, order] = sort(real(lambda), 'descend');
modes = struct();
modes.Lambda = lambda(order);
modes.Right = rightVectors(:, order);
modes.Left = leftVectors(:, order);
modes.ConditionNumber = conditionNumber(order);
modes.RightResidual = rightResidual(order);
modes.LeftResidual = leftResidual(order);
end

function [vectors, values] = local_dense_pairs(matrix, count, selector)
[allVectors, allValues] = eig(full(matrix), 'vector');
if ischar(selector) || isstring(selector)
    switch lower(char(selector))
        case {'largestreal','lr'}
            [~, order] = sort(real(allValues), 'descend');
        case {'smallestreal','sr'}
            [~, order] = sort(real(allValues), 'ascend');
        case {'largestabs','lm'}
            [~, order] = sort(abs(allValues), 'descend');
        otherwise
            error('Unsupported dense selector %s.', char(selector));
    end
else
    [~, order] = sort(abs(allValues - selector), 'ascend');
end
order = order(1:count);
vectors = allVectors(:, order);
values = diag(allValues(order));
end
