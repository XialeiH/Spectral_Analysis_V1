function out = GA_section5_from_response(response, varargin)
% Proposal section 5 covariance and reliability geometry from a response.
opts = struct('Dur', [], 'Dpr', [], 'Fano', 1, 'NoiseFloor', 1e-6, ...
    'Reg', 1e-6, 'EmpiricalTrials', [], 'Shrinkage', 0, ...
    'LowRankU', [], 'LowRankLambda', [], 'A', [], ...
    'CovarianceMode', 'linearized_discrete', 'SourceCovariance', [], ...
    'SolveCovariance', 'auto', 'MaxFullLyapN', 600, ...
    'NeumannMaxTerms', 100, 'NeumannTol', 1e-6);
opts = parse_opts(opts, varargin{:});

if ~isempty(opts.EmpiricalTrials)
    [Q, rhat] = GA_empirical_covariance(opts.EmpiricalTrials, opts.Shrinkage);
    out.SourceCovarianceModel = 'empirical_shrinkage';
    out.ResponseMean = rhat;
else
    r = GA_pack_state(response);
    if ~isempty(opts.LowRankU)
        diagVar = opts.Fano .* max(real(r), 0);
        Q = GA_lowrank_covariance(diagVar, opts.LowRankU, opts.LowRankLambda, opts.NoiseFloor);
        out.SourceCovarianceModel = 'lowrank_plus_diagonal';
    else
        Q = GA_rate_covariance(r, 'Fano', opts.Fano, 'NoiseFloor', opts.NoiseFloor);
        out.SourceCovarianceModel = 'diagonal_rate_dependent';
    end
    out.ResponseMean = r;
end

if ~isempty(opts.SourceCovariance)
    Q = opts.SourceCovariance;
    out.SourceCovarianceModel = 'user_supplied';
end

out.SourceCovariance = Q;
out.CovarianceMode = opts.CovarianceMode;
out.CovarianceEquation = '';
out.ResponseCovariance = [];
out.CovarianceSolved = false;
out.CovarianceConverged = false;
out.CovarianceSolution = '';
out.NeumannInfo = [];
out.GeometryCovarianceModel = '';
out.GeometryUsesSolvedCovariance = false;

switch lower(opts.CovarianceMode)
    case {'phase1', 'diagonal_rate_dependent'}
        out.CovarianceModel = out.SourceCovarianceModel;
        out.ResponseCovariance = Q;
        out.CovarianceSolved = true;
        out.CovarianceConverged = true;
    case {'linearized_discrete', 'phase2_discrete'}
        if isempty(opts.A)
            error('Phase-II discrete covariance requires option A, the fixed-point Jacobian.');
        end
        out.CovarianceModel = ['phase2_discrete_from_' out.SourceCovarianceModel];
        out.CovarianceEquation = 'Sigma = A*Sigma*A'' + Q';
        if should_solve_covariance(opts, size(opts.A, 1))
            try
                out.ResponseCovariance = GA_linearized_covariance(opts.A, Q, 'discrete');
                out.CovarianceSolved = true;
                out.CovarianceConverged = true;
                out.CovarianceSolution = 'direct_discrete_lyapunov';
            catch directErr
                out.CovarianceSolution = ['direct_failed: ' directErr.identifier];
            end
        end
        if isempty(out.ResponseCovariance)
            [out.ResponseCovariance, out.NeumannInfo] = GA_neumann_covariance(opts.A, Q, ...
                'MaxTerms', opts.NeumannMaxTerms, 'Tol', opts.NeumannTol);
            out.CovarianceSolved = true;
            out.CovarianceConverged = out.NeumannInfo.Converged && ~out.NeumannInfo.Diverged;
            out.CovarianceSolution = 'truncated_neumann_series';
        end
    case {'linearized_continuous', 'phase2_continuous'}
        if isempty(opts.A)
            error('Phase-II continuous covariance requires option A, the fixed-point Jacobian.');
        end
        out.CovarianceModel = ['phase2_continuous_from_' out.SourceCovarianceModel];
        out.CovarianceEquation = '(A-I)*Sigma + Sigma*(A-I)'' + Q = 0';
        if should_solve_covariance(opts, size(opts.A, 1))
            out.ResponseCovariance = GA_linearized_covariance(opts.A, Q, 'continuous');
            out.CovarianceSolved = true;
            out.CovarianceConverged = true;
            out.CovarianceSolution = 'direct_continuous_lyapunov';
        end
    otherwise
        error('Unknown CovarianceMode %s.', opts.CovarianceMode);
end

% Backward-compatible name; this is response Sigma only when solved.
out.RateCovariance = out.ResponseCovariance;

if ~isempty(opts.Dur)
    if out.CovarianceSolved
        SigmaGeom = out.ResponseCovariance;
        out.GeometryCovarianceModel = out.CovarianceModel;
        out.GeometryUsesSolvedCovariance = true;
    else
        SigmaGeom = out.SourceCovariance;
        out.GeometryCovarianceModel = out.SourceCovarianceModel;
    end
    out.Geometry = GA_reliability_geometry(opts.Dur, SigmaGeom, 'Dpr', opts.Dpr, 'Reg', opts.Reg);
else
    out.Geometry = [];
end
end

function tf = should_solve_covariance(opts, n)
if ischar(opts.SolveCovariance) || isstring(opts.SolveCovariance)
    mode = char(opts.SolveCovariance);
    switch lower(mode)
        case 'auto'
            tf = n <= opts.MaxFullLyapN;
        case {'true', 'on', 'yes'}
            tf = true;
        case {'false', 'off', 'no'}
            tf = false;
        otherwise
            error('Unknown SolveCovariance mode %s.', mode);
    end
else
    tf = logical(opts.SolveCovariance);
end
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
