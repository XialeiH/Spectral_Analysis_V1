function output = l6ns_phi_l6_ei(state,gain6,gainEI,context)
% Evaluate the moved-fixed-point map with L6 and I-to-E gains.
%
% gainEI scales only I->S and I->C. The I->I pathway is unchanged.

arguments
    state
    gain6 (1,1) double {mustBeFinite,mustBeNonnegative}
    gainEI (1,1) double {mustBeFinite,mustBeNonnegative}
    context (1,1) struct
end

state=state(:);
n=numel(state)/3;
if n~=round(n)
    error('L6EI:StateSize','State length must be divisible by three.');
end
s=state(1:n);
c=state(n+(1:n));
i=state(2*n+(1:n));
[sUse,cUse,iUse]=local_preprocess(s,c,i,context);

l4ES=(context.C_SS*sUse+context.C_SC*cUse)/context.L4SEp;
l4EC=(context.C_CS*sUse+context.C_CC*cUse)/context.L4CEp;
l4EI=(context.C_IS*sUse+context.C_IC*cUse)/context.L4IEp;

% beta_EI acts on inhibition received by the two excitatory populations.
l4IS=gainEI*(context.C_SI*iUse)/context.L4SIp;
l4IC=gainEI*(context.C_CI*iUse)/context.L4CIp;
l4II=(context.C_II*iUse)/context.L4IIp;

l6Input=gain6*l6ns_l6_dynamic(state,context,[1 1]);
outS=LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'S',l4ES,l4IS,context.ContrastUse,context.OrientationUse, ...
    context.PixLGNCtgr,l6Input);
outC=LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'C',l4EC,l4IC,context.ContrastUse,context.OrientationUse, ...
    context.PixLGNCtgr,l6Input);
outI=LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'I',l4EI,l4II,context.ContrastUse,context.OrientationUse, ...
    context.PixLGNCtgr,l6Input);
output=[outS(:);outC(:);outI(:)];
if isfield(context,'FixedPointCorrection') && ...
        ~isempty(context.FixedPointCorrection)
    output=output+context.FixedPointCorrection(:);
end
end

function [sUse,cUse,iUse]=local_preprocess(s,c,i,context)
if context.Isaturation
    eRaw=(1-context.CWeight)*s+context.CWeight*c;
    eBase=L6Convert(eRaw,context.EKpUse);
    adjustment=eBase./eRaw;
    adjustment(~isfinite(adjustment))=1;
    sUse=s.*adjustment;
    cUse=c.*adjustment;
    iUse=InhMulp(i,context.IKpUse);
else
    sUse=s;
    cUse=c;
    iUse=i;
end
end
