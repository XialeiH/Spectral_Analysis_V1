function fast = h96_nn_prepare_fast_config(ctx, contrast, orientation, modelCachePath)
% Precompute every condition-independent quantity used by the 50-step solve.
if nargin < 4 || isempty(modelCachePath)
    modelCachePath = fullfile(fileparts(mfilename('fullpath')), 'h96_models_double.mat');
end
cache = load(modelCachePath, 'models');

n = size(ctx.PixLGNCtgr, 1);
assert(size(ctx.PixLGNCtgr, 2) == 5, 'Expected five LGN categories.');

fast.n = n;
fast.fieldRows = ctx.N_HCOutY * ctx.NPixY;
fast.fieldCols = n / fast.fieldRows;
fast.kernel = ctx.L6Kernel;
fast.l6Max = 120;
fast.p = ctx.p;
fast.oneMinusP = 1 - ctx.p;
fast.EKp = ctx.EKpUse;
fast.IKp = ctx.IKpUse;
fast.ekpFast = h96_prepare_ekp_fast(ctx.EKpUse);
fast.ikpFast = h96_prepare_ikp_fast(ctx.IKpUse);
fast.ikpLookup = h96_prepare_ikp_lookup(ctx.IKpUse, 0.001);
fast.models = cache.models;
fast.PixLGNCtgr = ctx.PixLGNCtgr;

% Fuse nine recurrent sparse products into two sparse matrix-vector products.
fast.AE = [ctx.C_SS_meanU / ctx.L4SEp, ctx.C_SC_meanU / ctx.L4SEp; ...
           ctx.C_CS_meanU / ctx.L4CEp, ctx.C_CC_meanU / ctx.L4CEp; ...
           ctx.C_IS_mean  / ctx.L4IEp, ctx.C_IC_mean  / ctx.L4IEp];
fast.AI = [ctx.C_SI_mean / ctx.L4SIp; ...
           ctx.C_CI_mean / ctx.L4CIp; ...
           ctx.C_II_mean / ctx.L4IIp];

% Compile the current C1-smoothed L6 conversion to one reusable pp form.
l6Pars = ctx.L6pars{ctx.L6parId};
if iscell(l6Pars{end}) && numel(l6Pars{end}) >= 3 && strcmp(l6Pars{end}{1}, 'c1smooth')
    grid = l6Pars{end}{3};
    yGrid = L6Convert(grid, l6Pars{end}{2});
    fast.l6PP = pchip(grid, yGrid);
    fast.l6Pars = [];
else
    fast.l6PP = [];
    fast.l6Pars = l6Pars;
end

thetaPref = [0, 135, 90, 45, 0];
alpha = orientation(:);
if isscalar(alpha)
    alpha = repmat(alpha, n, 1);
end
assert(numel(alpha) == n, 'Orientation must be scalar or N-by-1.');

contrastVector = contrast(:);
if isscalar(contrastVector)
    contrastVector = repmat(contrastVector, n, 1);
end
assert(numel(contrastVector) == n, 'Contrast must be scalar or N-by-1.');

oriCos2 = zeros(n, 5);
contrastEff = repmat(contrastVector, 1, 5);
contrastEff(:,5) = 0;
for category = 1:5
    gamma = mod(alpha - thetaPref(category) + 90, 180) - 90;
    oriCos2(:,category) = cosd(2 .* abs(gamma));
end

% Category-major layout matches MATLAB's N-by-5 predictor output exactly.
fast.staticNormalized = struct();
for name = {'S','C','I'}
    population = name{1};
    model = fast.models.(population);
    x = zeros(n * 5, 5);
    x(:,1) = (oriCos2(:) - model.mu(1)) / model.sd(1);
    x(:,2) = (contrastEff(:) - model.mu(2)) / model.sd(2);
    fast.staticNormalized.(population) = x;
end
end
