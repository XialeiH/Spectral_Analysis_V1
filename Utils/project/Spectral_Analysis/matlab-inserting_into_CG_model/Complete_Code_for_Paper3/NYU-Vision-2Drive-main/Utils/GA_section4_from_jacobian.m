function out = GA_section4_from_jacobian(A, varargin)
% Proposal section 4 summary from the fixed-point update Jacobian A.
opts = struct('dPhi', [], 'FeatureNames', {{}}, 'ComputeFullResolvent', false);
opts = parse_opts(opts, varargin{:});

out.A = A;
out.L = A - speye(size(A, 1));
out.NonNormality = GA_nonnormality_index(A);

[V, D] = eig(full(A));
lambda = diag(D);
[~, idx] = sort(real(lambda), 'descend');
lambda = lambda(idx);
V = V(:, idx);

out.EigenValues = lambda;
out.RightEigenVectors = V;
out.MaxRealLambda = max(real(lambda));
out.StabilityMargin = 1 - out.MaxRealLambda;
out.SpectralRadius = max(abs(lambda));

if opts.ComputeFullResolvent
    n = size(A, 1);
    out.Resolvent = (eye(n) - full(A)) \ eye(n);
else
    out.Resolvent = [];
end

if ~isempty(opts.dPhi)
    out.TaskTangents = GA_resolvent_tangents(A, opts.dPhi);
    W = V' \ eye(size(V, 1));
    out.LeftEigenVectors = W;
    out.FeatureModeWeights = GA_feature_mode_weights(V, W, lambda, opts.dPhi);
    out.FeatureNames = opts.FeatureNames;
else
    out.TaskTangents = [];
    out.LeftEigenVectors = [];
    out.FeatureModeWeights = [];
    out.FeatureNames = {};
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
