function [LDEFixed, info] = GA_h96_response_map_fixedpoint(orientation, logSF, logTF, logContrast, ctx, varargin)
% Compute r(theta,logSF,logTF,logC) using the standard 100-step h96 iteration.
% The h96 local NN uses orientation and contrast directly. SF/TF are kept in
% info because the current h96 pipeline fixes the LGN drive upstream.
if nargin < 5 || isempty(ctx)
    ctx = GA_h96_context_from_base();
end

opts = struct('InitialState', [], 'IterationEpoch', 100);
opts = parse_opts(opts, varargin{:});

contrast = exp(logContrast);
sf = exp(logSF);
tf = exp(logTF);

if isempty(opts.InitialState)
    iniState = h96_initial_state(ctx);
    initialStateSource = 'standard_h96_initial_state';
else
    initialStateSource = 'user_supplied';
    iniState = opts.InitialState;
end

contrastUse = contrast * ones(size(iniState.S));
orientationUse = orientation * ones(size(iniState.S));

iterTimer = tic;
[LDEItr,~,~,~,~,~,NANFlag] = ...
    LDEIteration_135FuncMain_CombDom_RealLGNL6_MLP6D_prefAngle(...
    ctx.PixLGNCtgr, ctx.L6Kernel, iniState, ctx.p, ctx.L6pars{ctx.L6parId}, opts.IterationEpoch,...
    ctx.C_SS_meanU, ctx.C_CS_meanU, ctx.C_IS_mean,...
    ctx.C_SC_meanU, ctx.C_CC_meanU, ctx.C_IC_mean,...
    ctx.C_SI_mean, ctx.C_CI_mean, ctx.C_II_mean,...
    ctx.L4SEp, ctx.L4SIp, ...
    ctx.L4CEp, ctx.L4CIp, ...
    ctx.L4IEp, ctx.L4IIp, ...
    ctx.L4EmeshXAll, ctx.L4ImeshYAll, ctx.LDEFrfuncAll, ...
    contrastUse, orientationUse, ...
    ctx.N_HCOutY, ctx.NPixX, ctx.NPixY, ctx.Isaturation, 'xn', ctx.EKpUse, ctx.IKpUse);
iterationTime = toc(iterTimer);

LDEFixed = LDEItr{end};
xFixed = GA_pack_state(LDEFixed);

info = struct();
info.orientation = orientation;
info.logSF = logSF;
info.logTF = logTF;
info.logContrast = logContrast;
info.SF = sf;
info.TF = tf;
info.Contrast = contrast;
info.InitialStateSource = initialStateSource;
info.IterationEpoch = opts.IterationEpoch;
info.IterationTime = iterationTime;
info.NANFlag = NANFlag;
info.SF_TF_note = 'Current h96 baseline fixes SF/TF upstream in the LGN drive; this solve uses the current ctx LGN map.';
info.ResponseVector = xFixed;
end

function state = h96_initial_state(ctx)
if isfield(ctx, 'IniTest')
    state = ctx.IniTest;
    return
end
n = size(ctx.PixLGNCtgr, 1);
state = struct();
state.S = 2.5 * ones(n, 1);
state.C = 8.0 * ones(n, 1);
state.I = 18.0 * ones(n, 1);
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
