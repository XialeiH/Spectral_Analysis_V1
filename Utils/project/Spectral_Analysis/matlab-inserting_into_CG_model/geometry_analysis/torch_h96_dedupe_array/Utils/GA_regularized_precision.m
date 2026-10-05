function Prec = GA_regularized_precision(Sigma, reg)
% Stable precision matrix for geometry calculations.
if nargin < 2 || isempty(reg)
    reg = 1e-6;
end
n = size(Sigma, 1);
if issparse(Sigma)
    [ii, jj, vv] = find(Sigma);
    if all(ii == jj)
        if isempty(ii)
            d = zeros(n, 1);
        else
            d = accumarray(ii, vv, [n 1], @sum, 0);
        end
        Prec = spdiags(1 ./ (d + reg), 0, n, n);
        return
    end
end
SigmaUse = full(0.5 * (Sigma + Sigma')) + reg * eye(n);
Prec = SigmaUse \ eye(n);
Prec = 0.5 * (Prec + Prec');
end
