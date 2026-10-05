function geom = GA_reliability_geometry(Dur, Sigma, varargin)
% Compute g, H, C blocks from proposal section 5.9 / equation 71.
opts = struct('Dpr', [], 'Reg', 1e-6);
opts = parse_opts(opts, varargin{:});

Prec = GA_regularized_precision(Sigma, opts.Reg);
geom.Precision = Prec;
geom.g = Dur' * Prec * Dur;

if ~isempty(opts.Dpr)
    geom.H = opts.Dpr' * Prec * opts.Dpr;
    geom.C = Dur' * Prec * opts.Dpr;
else
    geom.H = [];
    geom.C = [];
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
