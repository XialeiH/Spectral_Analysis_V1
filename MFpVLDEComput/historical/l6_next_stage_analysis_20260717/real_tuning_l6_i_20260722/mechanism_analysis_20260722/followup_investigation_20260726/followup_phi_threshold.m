function output = followup_phi_threshold(state,context,gain6,gainI,mode, ...
        edgeScale,blendScale)
% True-tuning map with raw or threshold-shifted extension response.

if mode=="raw"
    output=mechanism_phi_variant(state,context,gain6,gainI,'raw');
    return
end
inputs=mechanism_state_inputs(state,context,gain6,gainI);
outS=local_aggregate('S',inputs.S.L4E,inputs.S.L4I,inputs.L6,context, ...
    edgeScale,blendScale);
outC=local_aggregate('C',inputs.C.L4E,inputs.C.L4I,inputs.L6,context, ...
    edgeScale,blendScale);
outI=local_aggregate('I',inputs.I.L4E,inputs.I.L4I,inputs.L6,context, ...
    edgeScale,blendScale);
output=[outS(:);outC(:);outI(:)];
if isfield(context,'FixedPointCorrection') && ~isempty(context.FixedPointCorrection)
    output=output+context.FixedPointCorrection(:);
end
end

function output=local_aggregate(celltype,E,I,L6,context,edgeScale,blendScale)
output=zeros(size(E));
thetaMap=[0 135 90 45];
for lgn=1:5
    contrast=context.ContrastUse;
    oriCos2=ones(size(context.OrientationUse));
    if lgn<=4
        gamma=abs(mod(double(context.OrientationUse)-thetaMap(lgn)+90,180)-90);
        oriCos2=cosd(2*gamma);
    else
        contrast(:)=0;
    end
    y=followup_domain_extended_response(celltype,E,I,L6,contrast,oriCos2, ...
        edgeScale,blendScale);
    output=output+context.PixLGNCtgr(:,lgn).*y;
end
end
