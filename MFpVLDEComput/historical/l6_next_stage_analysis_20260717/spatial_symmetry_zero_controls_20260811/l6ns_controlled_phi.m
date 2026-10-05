function output=l6ns_controlled_phi(state,context,operators)
% One exact h96 response-map evaluation using controlled spatial operators.

n=prod(context.MapSize);
state=state(:);
s=state(1:n); c=state(n+(1:n)); i=state(2*n+(1:n));
if context.Isaturation
    eRaw=(1-context.CWeight)*s+context.CWeight*c;
    adjustment=L6Convert(eRaw,context.EKpUse)./eRaw;
    adjustment(~isfinite(adjustment))=1;
    sUse=s.*adjustment;
    cUse=c.*adjustment;
    iUse=InhMulp(i,context.IKpUse);
else
    sUse=s; cUse=c; iUse=i;
end

l4ES=(operators.C_SS*sUse+operators.C_SC*cUse)/context.L4SEp;
l4EC=(operators.C_CS*sUse+operators.C_CC*cUse)/context.L4CEp;
l4EI=(operators.C_IS*sUse+operators.C_IC*cUse)/context.L4IEp;
l4IS=(operators.C_SI*iUse)/context.L4SIp;
l4IC=(operators.C_CI*iUse)/context.L4CIp;
l4II=(operators.C_II*iUse)/context.L4IIp;

eUse=(1-context.CWeight)*sUse+context.CWeight*cUse;
l6=L6Convert(reshape(operators.L6*eUse,context.MapSize), ...
    context.L6Parameters);
l6=min(max(l6(:)/3,1),40);

sOut=LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'S',l4ES,l4IS,context.ContrastUse,context.OrientationUse, ...
    context.PixLGNCtgr,l6);
cOut=LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'C',l4EC,l4IC,context.ContrastUse,context.OrientationUse, ...
    context.PixLGNCtgr,l6);
iOut=LDEIterFunc_Grating_135Func_RealLGNL6_MLP6D_prefAngle( ...
    'I',l4EI,l4II,context.ContrastUse,context.OrientationUse, ...
    context.PixLGNCtgr,l6);
output=[sOut(:);cOut(:);iOut(:)];
end
