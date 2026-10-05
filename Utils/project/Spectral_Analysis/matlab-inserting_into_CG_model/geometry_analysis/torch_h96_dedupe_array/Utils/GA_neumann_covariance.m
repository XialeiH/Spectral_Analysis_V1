function [Sigma, info] = GA_neumann_covariance(A, Q, varargin)
% Truncated Neumann-series covariance for Sigma = A*Sigma*A' + Q.
opts = struct('MaxTerms', 100, 'Tol', 1e-6, 'DivergenceFactor', 1e12);
opts = parse_opts(opts, varargin{:});

Sigma = Q;
Term = Q;
baseNorm = max(norm(Q, 'fro'), eps);
sigmaNorm = max(norm(Sigma, 'fro'), eps);

info = struct();
info.Method = 'truncated_neumann_series';
info.Equation = 'Sigma = sum_{k=0}^{K} A^k Q (A'')^k';
info.MaxTerms = opts.MaxTerms;
info.Tol = opts.Tol;
info.Converged = false;
info.TermsUsed = 0;
info.LastTermNorm = norm(Term, 'fro');
info.LastRelativeTermNorm = info.LastTermNorm / sigmaNorm;
info.Diverged = false;

for k = 1:opts.MaxTerms
    Term = A * Term * A';
    Term = 0.5 * (Term + Term');
    termNorm = norm(Term, 'fro');

    if ~isfinite(termNorm) || termNorm > opts.DivergenceFactor * baseNorm
        info.Diverged = true;
        info.TermsUsed = k - 1;
        break
    end

    Sigma = Sigma + Term;
    sigmaNorm = max(norm(Sigma, 'fro'), eps);
    relTermNorm = termNorm / sigmaNorm;

    info.TermsUsed = k;
    info.LastTermNorm = termNorm;
    info.LastRelativeTermNorm = relTermNorm;

    if relTermNorm < opts.Tol
        info.Converged = true;
        break
    end
end

Sigma = 0.5 * (Sigma + Sigma');
end

function opts = parse_opts(opts, varargin)
if mod(numel(varargin), 2) ~= 0
    error('Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    value = varargin{k+1};
    if ~isfield(opts, name)
        error('Unknown option %s.', name);
    end
    opts.(name) = value;
end
end
