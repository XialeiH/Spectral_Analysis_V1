function out = GA_fixed_point_summary(response, varargin)
% Proposal section 4.1 fixed-point response summary.
opts = struct('NextResponse', [], 'NormWeights', []);
opts = parse_opts(opts, varargin{:});

out.Response = response;
out.ResponseVector = GA_pack_state(response);
out.ResponseNorm = norm(out.ResponseVector);

if ~isempty(opts.NextResponse)
    nextVec = GA_pack_state(opts.NextResponse);
    out.ResidualVector = nextVec - out.ResponseVector;
    out.ResidualNorm = norm(out.ResidualVector);
else
    out.ResidualVector = [];
    out.ResidualNorm = [];
end

if ~isempty(opts.NextResponse) && ~isempty(opts.NormWeights) && isstruct(response)
    w = opts.NormWeights;
    out.HCResidualNorm = HC_norm_diff(response.S, response.C, response.I, ...
        opts.NextResponse.S, opts.NextResponse.C, opts.NextResponse.I, ...
        w.wC, w.aE, w.aI);
else
    out.HCResidualNorm = [];
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
