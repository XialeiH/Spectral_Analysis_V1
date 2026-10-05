function Sigma = GA_rate_covariance(r, varargin)
% Diagonal rate-dependent covariance from proposal section 5.5.
opts = struct('Fano', 1, 'NoiseFloor', 1e-6);
opts = parse_opts(opts, varargin{:});

r = max(real(r(:)), 0);
if isscalar(opts.Fano)
    fano = opts.Fano * ones(size(r));
else
    fano = opts.Fano(:);
end
Sigma = spdiags(fano .* r + opts.NoiseFloor, 0, numel(r), numel(r));
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
