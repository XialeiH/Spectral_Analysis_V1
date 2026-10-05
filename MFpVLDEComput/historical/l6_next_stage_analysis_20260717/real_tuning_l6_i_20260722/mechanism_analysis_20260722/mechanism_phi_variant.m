function output = mechanism_phi_variant(state,context,gain6,gainI,responseMode)
% Evaluate true tuning with either the exact extension or the raw h96 NN.

inputs = mechanism_state_inputs(state,context,gain6,gainI);
outS = mechanism_aggregate_response('S',inputs.S.L4E,inputs.S.L4I, ...
    inputs.L6,context,responseMode);
outC = mechanism_aggregate_response('C',inputs.C.L4E,inputs.C.L4I, ...
    inputs.L6,context,responseMode);
outI = mechanism_aggregate_response('I',inputs.I.L4E,inputs.I.L4I, ...
    inputs.L6,context,responseMode);
output = [outS(:);outC(:);outI(:)];
if isfield(context,'FixedPointCorrection') && ~isempty(context.FixedPointCorrection)
    output = output+context.FixedPointCorrection(:);
end
end
