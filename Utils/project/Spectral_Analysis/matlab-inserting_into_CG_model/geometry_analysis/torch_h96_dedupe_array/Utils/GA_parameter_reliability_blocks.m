function out = GA_parameter_reliability_blocks(g, H, C, paramNames)
% Derived P-geometry blocks from proposal sections 8.1, 8.2, and 13.3.
out = struct();
out.g = g;
out.H = H;
out.C = C;
out.ParameterNames = paramNames;

if isempty(g) || isempty(H) || isempty(C)
    out.B = [];
    out.Hperp = [];
    out.NormalizedConfound = [];
    return;
end

out.B = safe_solve(g, C);
out.Hperp = H - C' * safe_solve(g, C);

nTask = size(C, 1);
nParam = size(C, 2);
rho = nan(nTask, nParam);
for i = 1:nTask
    gi = real(g(i, i));
    for a = 1:nParam
        ha = real(H(a, a));
        denom = sqrt(max(gi, 0) * max(ha, 0));
        if denom > 0
            rho(i, a) = abs(C(i, a)) / denom;
        end
    end
end
out.NormalizedConfound = rho;
out.HDiagonal = real(diag(H));
out.HperpDiagonal = real(diag(out.Hperp));
out.ResidualFraction = out.HperpDiagonal ./ out.HDiagonal;
out.ResidualFraction(~isfinite(out.ResidualFraction)) = NaN;
end

function X = safe_solve(A, B)
try
    X = A \ B;
catch
    X = pinv(full(A)) * B;
end
end
