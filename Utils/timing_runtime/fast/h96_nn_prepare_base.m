function base = h96_nn_prepare_base(ctx, precision, modelCachePath)
% Prepare condition-independent h96 NN fixed-point data once.
if nargin<2 || isempty(precision)
    precision = 'double';
end
if nargin<3 || isempty(modelCachePath)
    modelCachePath = fullfile(fileparts(mfilename('fullpath')), 'h96_models_double.mat');
end
assert(strcmp(precision, 'double') || strcmp(precision, 'single'), ...
    'Precision must be double or single.');

cache = load(modelCachePath, 'models');
n = size(ctx.PixLGNCtgr, 1);
assert(size(ctx.PixLGNCtgr, 2)==5, 'Expected five LGN categories.');

base.n = n;
base.fieldRows = ctx.N_HCOutY*ctx.NPixY;
base.fieldCols = n/base.fieldRows;
base.kernel = ctx.L6Kernel;
base.l6Max = 120;
base.p = ctx.p;
base.oneMinusP = 1-ctx.p;
base.EKp = ctx.EKpUse;
base.IKp = ctx.IKpUse;
base.nnPrecision = precision;

% The recurrent operators are invariant across contrast/orientation conditions.
base.AE = [ctx.C_SS_meanU/ctx.L4SEp, ctx.C_SC_meanU/ctx.L4SEp; ...
           ctx.C_CS_meanU/ctx.L4CEp, ctx.C_CC_meanU/ctx.L4CEp; ...
           ctx.C_IS_mean /ctx.L4IEp, ctx.C_IC_mean /ctx.L4IEp];
base.AI = [ctx.C_SI_mean/ctx.L4SIp; ...
           ctx.C_CI_mean/ctx.L4CIp; ...
           ctx.C_II_mean/ctx.L4IIp];

l6Pars = ctx.L6pars{ctx.L6parId};
if iscell(l6Pars{end}) && numel(l6Pars{end})>=3 && strcmp(l6Pars{end}{1}, 'c1smooth')
    grid = l6Pars{end}{3};
    base.l6PP = pchip(grid, L6Convert(grid, l6Pars{end}{2}));
    base.l6Pars = [];
else
    base.l6PP = [];
    base.l6Pars = l6Pars;
end

[rowIndex, categoryIndex, pairWeight] = find(ctx.PixLGNCtgr);
pairCount = numel(rowIndex);
base.mixRows = rowIndex;
base.mixCategories = categoryIndex;
base.mixAggregate = sparse(rowIndex, 1:pairCount, pairWeight, n, pairCount);
base.activePairCount = pairCount;
base.thetaPref = [0, 135, 90, 45, 0];

names = {'S', 'C', 'I'};
base.mu = zeros(1, 5, 3, precision);
base.sd = zeros(1, 5, 3, precision);
base.W1T = zeros(5, 96, 3, precision);
base.b1 = zeros(1, 96, 3, precision);
base.W2T = zeros(96, 96, 3, precision);
base.b2 = zeros(1, 96, 3, precision);
base.W3T = zeros(96, 96, 3, precision);
base.b3 = zeros(1, 96, 3, precision);
base.W4T = zeros(96, 1, 3, precision);
base.b4 = zeros(1, 1, 3, precision);
for population = 1:3
    model = cache.models.(names{population});
    base.mu(:,:,population) = cast(model.mu, precision);
    base.sd(:,:,population) = cast(model.sd, precision);
    base.W1T(:,:,population) = cast(model.W1.', precision);
    base.b1(:,:,population) = cast(model.b1, precision);
    base.W2T(:,:,population) = cast(model.W2.', precision);
    base.b2(:,:,population) = cast(model.b2, precision);
    base.W3T(:,:,population) = cast(model.W3.', precision);
    base.b3(:,:,population) = cast(model.b3, precision);
    base.W4T(:,:,population) = cast(model.W4.', precision);
    base.b4(:,:,population) = cast(model.b4, precision);
end
end
