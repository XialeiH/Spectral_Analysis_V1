function output = l6ns_phi(state, l6Weight, context, sourceScale, inhibitoryGain, interventionMode)
% Evaluate Phi with fixed-point-preserving L6 and inhibitory pathway gains.
%
% l6Weight may be scalar or [wS wC wI].  sourceScale is [betaS betaC]
% and acts only inside the L6 feedback path, around the stored fixed point.

if nargin < 4 || isempty(sourceScale)
    sourceScale = [1 1];
end
if nargin < 5 || isempty(inhibitoryGain)
    inhibitoryGain = 1;
end
if nargin < 6 || isempty(interventionMode)
    interventionMode = 'fpp';
end
if ~any(strcmp(interventionMode,{'fpp','true'}))
    error('interventionMode must be fpp or true.');
end
if ~isscalar(inhibitoryGain) || ~isfinite(inhibitoryGain)
    error('inhibitoryGain must be one finite scalar.');
end
if isstruct(state)
    s = state.S(:); c = state.C(:); i = state.I(:);
else
    state = state(:);
    n = numel(state) / 3;
    s = state(1:n); c = state(n+(1:n)); i = state(2*n+(1:n));
end
if isscalar(l6Weight)
    l6Weight = repmat(l6Weight, 1, 3);
end
if numel(l6Weight) ~= 3
    error('l6Weight must be scalar or [wS wC wI].');
end
if numel(sourceScale) ~= 2
    error('sourceScale must be [betaS betaC].');
end

[sUse,cUse,iUse] = local_preprocess(s,c,i,context);
l4ES = (context.C_SS * sUse + context.C_SC * cUse) / context.L4SEp;
l4EC = (context.C_CS * sUse + context.C_CC * cUse) / context.L4CEp;
l4EI = (context.C_IS * sUse + context.C_IC * cUse) / context.L4IEp;
l4IS = (context.C_SI * iUse) / context.L4SIp;
l4IC = (context.C_CI * iUse) / context.L4CIp;
l4II = (context.C_II * iUse) / context.L4IIp;

% Preserve the inhibitory pathway value at f* while multiplying its
% complete I-source derivative by inhibitoryGain.
if inhibitoryGain ~= 1 && strcmp(interventionMode,'fpp')
    fixed = context.FixedPoint(:);
    n = numel(fixed) / 3;
    [~,~,iFixedUse] = local_preprocess(fixed(1:n), ...
        fixed(n+(1:n)),fixed(2*n+(1:n)),context);
    l4ISFixed = (context.C_SI * iFixedUse) / context.L4SIp;
    l4ICFixed = (context.C_CI * iFixedUse) / context.L4CIp;
    l4IIFixed = (context.C_II * iFixedUse) / context.L4IIp;
    l4IS = l4ISFixed + inhibitoryGain * (l4IS-l4ISFixed);
    l4IC = l4ICFixed + inhibitoryGain * (l4IC-l4ICFixed);
    l4II = l4IIFixed + inhibitoryGain * (l4II-l4IIFixed);
elseif inhibitoryGain ~= 1
    l4IS = inhibitoryGain*l4IS;
    l4IC = inhibitoryGain*l4IC;
    l4II = inhibitoryGain*l4II;
end

% Compensated source ablation: only the L6 pathway sees the scaled
% displacement from the fixed point, so Phi(f*) is unchanged.
l6Dynamic = l6ns_l6_dynamic([s;c;i],context,sourceScale);

if strcmp(interventionMode,'true')
    l6S = (1-l6Weight(1))*l6Dynamic;
    l6C = (1-l6Weight(2))*l6Dynamic;
    l6I = (1-l6Weight(3))*l6Dynamic;
else
    l6S = l6Weight(1) * context.FixedL6LibInd + (1-l6Weight(1)) * l6Dynamic;
    l6C = l6Weight(2) * context.FixedL6LibInd + (1-l6Weight(2)) * l6Dynamic;
    l6I = l6Weight(3) * context.FixedL6LibInd + (1-l6Weight(3)) * l6Dynamic;
end

outS = LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'S',l4ES,l4IS,context.ContrastUse,context.OrientationUse,context.PixLGNCtgr,l6S);
outC = LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'C',l4EC,l4IC,context.ContrastUse,context.OrientationUse,context.PixLGNCtgr,l6C);
outI = LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'I',l4EI,l4II,context.ContrastUse,context.OrientationUse,context.PixLGNCtgr,l6I);
output = [outS(:); outC(:); outI(:)];
if isfield(context,'FixedPointCorrection') && ~isempty(context.FixedPointCorrection)
    output = output + context.FixedPointCorrection(:);
end
end

function [sUse,cUse,iUse] = local_preprocess(s,c,i,context)
if context.Isaturation
    eRaw = (1-context.CWeight) * s + context.CWeight * c;
    eBase = L6Convert(eRaw, context.EKpUse);
    adjustment = eBase ./ eRaw;
    adjustment(~isfinite(adjustment)) = 1;
    sUse = s .* adjustment;
    cUse = c .* adjustment;
    iUse = InhMulp(i, context.IKpUse);
else
    sUse = s;
    cUse = c;
    iUse = i;
end
end
