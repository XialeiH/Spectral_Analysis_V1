function inputs = mechanism_state_inputs(state, context, gain6, gainI)
% Reconstruct the exact local-response inputs used by true pathway tuning.

state = state(:);
n = numel(state)/3;
s = state(1:n);
c = state(n+(1:n));
i = state(2*n+(1:n));
wC = context.CWeight;

if context.Isaturation
    eRaw = (1-wC)*s+wC*c;
    eBase = L6Convert(eRaw,context.EKpUse);
    adjustment = eBase./eRaw;
    adjustment(~isfinite(adjustment)) = 1;
    sUse = s.*adjustment;
    cUse = c.*adjustment;
    iUse = InhMulp(i,context.IKpUse);
else
    sUse = s;
    cUse = c;
    iUse = i;
end

inputs.S.L4E = (context.C_SS*sUse+context.C_SC*cUse)/context.L4SEp;
inputs.C.L4E = (context.C_CS*sUse+context.C_CC*cUse)/context.L4CEp;
inputs.I.L4E = (context.C_IS*sUse+context.C_IC*cUse)/context.L4IEp;
inputs.S.L4I = gainI*(context.C_SI*iUse)/context.L4SIp;
inputs.C.L4I = gainI*(context.C_CI*iUse)/context.L4CIp;
inputs.I.L4I = gainI*(context.C_II*iUse)/context.L4IIp;
inputs.L6Base = l6ns_l6_dynamic(state,context,[1 1]);
inputs.L6 = gain6*inputs.L6Base;
inputs.SUse = sUse;
inputs.CUse = cUse;
inputs.IUse = iUse;
end
