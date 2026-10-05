function correction = h96_pref6D_strength_moment_correction( ...
    celltype, L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse, PixLGNCtgr, Perturb)
% Second-cumulant correction for the no-new-library S_II approximation.

correction = zeros(size(L4IUse));
if ~getfield_default(Perturb, 'MomentCorrection', false)
    return;
end

momentParam = getfield_default(Perturb, 'MomentParam', '');
if ~strcmp(momentParam, 'S_II')
    error('h96_pref6D_strength_moment_correction:UnsupportedParameter', ...
        'Moment correction is currently implemented only for S_II, not %s.', momentParam);
end
if ~strcmp(celltype, 'I')
    return;
end

scale = getfield_default(Perturb, 'S_II', 1);
if scale == 1
    return;
end

% L4I is measured in events/s. Reproduce the normalized, truncated GABA
% kernel used by LIFo_SinglePixel_CurrDecomp and convert its H2 from ms to
% the NN's Hz coordinate: eta = H2 / (1e-3 events/ms per events/s).
dt = 0.1;
tauR = 0.5;
tauD = 5;
t = 0:dt:(3 * tauD);
kernel = (exp(-t / tauD) - exp(-t / tauR)) / (tauD - tauR);
kernel = kernel / (sum(kernel) * dt);
etaHz = (sum(kernel.^2) * dt) / 1e-3;

% Estimate d2 Phi_I / dL4I2 from the existing analytic first derivative.
% The default step is 0.005 of the h96 L4I training standard deviation.
curvatureStep = getfield_default(Perturb, 'MomentCurvatureStep', 0.005 * 7195.099609);
if ~isscalar(curvatureStep) || ~isfinite(curvatureStep) || curvatureStep <= 0
    error('h96_pref6D_strength_moment_correction:InvalidStep', ...
        'MomentCurvatureStep must be a positive finite scalar.');
end

curvature = centered_curvature_from_grad(celltype, L4EUse, L4IUse, L6Use, ...
    ContrastUse, OrientationUse, PixLGNCtgr, curvatureStep);
if getfield_default(Perturb, 'MomentRichardson', true)
    curvatureHalf = centered_curvature_from_grad(celltype, L4EUse, L4IUse, L6Use, ...
        ContrastUse, OrientationUse, PixLGNCtgr, curvatureStep / 2);
    curvature = (4 * curvatureHalf - curvature) / 3;
end

% L4IUse is already the mean-matched coordinate x' = scale*x because the
% connectivity matrix has been scaled. Therefore
% (scale^2-scale)*x = (scale-1)*x'.
correction = 0.5 * etaHz * (scale - 1) .* L4IUse .* curvature;
correction(~isfinite(correction)) = 0;
end

function curvature = centered_curvature_from_grad( ...
    celltype, L4EUse, L4IUse, L6Use, ContrastUse, OrientationUse, PixLGNCtgr, step)
[~, gradPlus, ~] = local_h96baseline_pref6D_grads( ...
    celltype, L4EUse, L4IUse + step, L6Use, ...
    ContrastUse, OrientationUse, PixLGNCtgr);
[~, gradMinus, ~] = local_h96baseline_pref6D_grads( ...
    celltype, L4EUse, L4IUse - step, L6Use, ...
    ContrastUse, OrientationUse, PixLGNCtgr);
curvature = (gradPlus - gradMinus) / (2 * step);
end

function val = getfield_default(s, name, defaultVal)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end
