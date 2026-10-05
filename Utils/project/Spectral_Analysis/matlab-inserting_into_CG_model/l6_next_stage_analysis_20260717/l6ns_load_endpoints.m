function data = l6ns_load_endpoints(cfg)
% Load the dynamic-L6 and fixed-L6 endpoint Jacobians for one condition.

fileW0 = local_weight_file(cfg.GeometryRoot, cfg.Angle, cfg.Contrast, 0);
fileW1 = local_weight_file(cfg.GeometryRoot, cfg.Angle, cfg.Contrast, 1);

s0 = load(fileW0, 'Section4', 'GeometryMetrics');
s1 = load(fileW1, 'Section4', 'GeometryMetrics');

j0 = local_get_jacobian(s0.Section4, fileW0);
a = local_get_jacobian(s1.Section4, fileW1);
b = j0 - a;

data = struct();
data.A = sparse(a);
data.B = sparse(b);
data.J0 = sparse(j0);
data.FileW0 = fileW0;
data.FileW1 = fileW1;
data.Angle = cfg.Angle;
data.Contrast = cfg.Contrast;
data.Dimension = size(a, 1);
if mod(data.Dimension, 3) ~= 0
    error('Expected a three-population Jacobian, got dimension %d.', data.Dimension);
end
data.PopulationSize = data.Dimension / 3;
data.MapSize = local_map_size(s0, data.PopulationSize);
data.ExcitatoryCWeight = cfg.ExcitatoryCWeight;
data.FixedPoint = local_fixed_point(s0);
data.FixedPointEndpointDifference = local_fixed_point_difference(s0, s1);
data.EndpointAffineError = norm((data.A + data.B) - data.J0, 'fro') / max(norm(data.J0, 'fro'), eps);
data.L6Kernel = local_l6_kernel(s0);

fprintf('Loaded A=J(1) and B=J(0)-J(1) for angle %.2f contrast %g.\n', ...
    data.Angle, data.Contrast);
fprintf('n=%d, nnz(A)=%d, nnz(B)=%d, fixed-point endpoint difference %.3e.\n', ...
    data.Dimension, nnz(data.A), nnz(data.B), data.FixedPointEndpointDifference);
end

function kernel = local_l6_kernel(s)
kernel = [];
if isfield(s, 'GeometryMetrics') && isfield(s.GeometryMetrics, 'JacobianInputs') ...
        && isfield(s.GeometryMetrics.JacobianInputs, 'L6Kernel')
    kernel = s.GeometryMetrics.JacobianInputs.L6Kernel;
end
end

function fileName = local_weight_file(rootDir, angleValue, contrastValue, weight)
raw = sprintf('L6eqW%.2f', weight);
tags = {strrep(raw, '.', 'p'), strrep(strrep(raw, '.', 'p'), '-', 'm')};
for i = 1:numel(tags)
    candidate = fullfile(rootDir, sprintf( ...
        'geometry_sections4_5_h96baseline_%s_contr%d_angle_%.2f.mat', ...
        tags{i}, contrastValue, angleValue));
    if isfile(candidate)
        fileName = candidate;
        return
    end
end
error('Could not find the L6-weight file for w=%.6g under %s.', weight, rootDir);
end

function j = local_get_jacobian(section4, fileName)
if isfield(section4, 'A')
    j = section4.A;
elseif isfield(section4, 'Jacobian')
    j = section4.Jacobian;
else
    error('No Jacobian field was found in %s.', fileName);
end
end

function mapSize = local_map_size(s, populationSize)
if isfield(s, 'GeometryMetrics') && isfield(s.GeometryMetrics, 'FixedPointMaps') ...
        && isfield(s.GeometryMetrics.FixedPointMaps, 'S')
    mapSize = size(s.GeometryMetrics.FixedPointMaps.S);
else
    side = round(sqrt(populationSize));
    if side * side ~= populationSize
        error('Cannot infer a two-dimensional map size from %d pixels.', populationSize);
    end
    mapSize = [side side];
end
end

function fixedPoint = local_fixed_point(s)
fixedPoint = [];
if isfield(s, 'LDEfixedpoint')
    fixedPoint = [s.LDEfixedpoint.S(:); s.LDEfixedpoint.C(:); s.LDEfixedpoint.I(:)];
elseif isfield(s, 'GeometryMetrics') && isfield(s.GeometryMetrics, 'FixedPoint') ...
        && isfield(s.GeometryMetrics.FixedPoint, 'ResponseVector')
    fixedPoint = s.GeometryMetrics.FixedPoint.ResponseVector(:);
elseif isfield(s, 'GeometryMetrics') && isfield(s.GeometryMetrics, 'FixedPointMaps')
    fp = s.GeometryMetrics.FixedPointMaps;
    fixedPoint = [fp.S(:); fp.C(:); fp.I(:)];
end
end

function difference = local_fixed_point_difference(s0, s1)
f0 = local_fixed_point(s0);
f1 = local_fixed_point(s1);
if isempty(f0) || isempty(f1) || numel(f0) ~= numel(f1)
    difference = NaN;
else
    difference = norm(f0 - f1) / max(norm(f0), eps);
end
end
