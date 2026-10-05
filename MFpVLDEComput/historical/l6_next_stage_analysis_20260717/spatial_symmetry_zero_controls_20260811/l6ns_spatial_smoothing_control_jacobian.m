function J = l6ns_spatial_smoothing_control_jacobian(setup, l4Alpha, l6Alpha)
% Replace spatial smoothing by matched-row-sum local coupling at fixed gains.

context = setup.Context;
state = context.FixedPoint(:);
n = numel(state)/3;
s = state(1:n);
c = state(n+(1:n));
i = state(2*n+(1:n));

[sUse,cUse,iUse,dSuseS,dSuseC,dCuseS,dCuseC,dIuseI] = ...
    local_preprocess_with_derivative(s,c,i,context);

% Evaluate all local NN and saturation derivatives on the canonical model.
l4ES = (context.C_SS*sUse + context.C_SC*cUse)/context.L4SEp;
l4EC = (context.C_CS*sUse + context.C_CC*cUse)/context.L4CEp;
l4EI = (context.C_IS*sUse + context.C_IC*cUse)/context.L4IEp;
l4IS = (context.C_SI*iUse)/context.L4SIp;
l4IC = (context.C_CI*iUse)/context.L4CIp;
l4II = (context.C_II*iUse)/context.L4IIp;

wS = 1-context.CWeight;
wC = context.CWeight;
eUse = wS*sUse+wC*cUse;
kBase = build_conv_matrix_circular(context.L6Kernel, ...
    context.MapSize(1),context.MapSize(2));
convolved = kBase*eUse;
l6Raw = L6Convert(reshape(convolved,context.MapSize),context.L6Parameters);
l6Raw = l6Raw(:)/3;
l6Used = min(max(l6Raw,1),40);
l6Gradient = L6Convert_grad(reshape(convolved,context.MapSize), ...
    context.L6Parameters);
l6Gradient = l6Gradient(:);
l6Mask = l6Raw>1 & l6Raw<40;
dL6Base = spdiags((l6Mask/3).*l6Gradient,0,n,n);

[dS4E,dS4I,dS6] = local_h96baseline_pref6D_grads('S',l4ES,l4IS,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
[dC4E,dC4I,dC6] = local_h96baseline_pref6D_grads('C',l4EC,l4IC,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
[dI4E,dI4I,dI6] = local_h96baseline_pref6D_grads('I',l4EI,l4II,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);

% Change only the spatial propagation operators in D Phi.
cSS = local_flatten_smoothing(context.C_SS,l4Alpha);
cSC = local_flatten_smoothing(context.C_SC,l4Alpha);
cSI = local_flatten_smoothing(context.C_SI,l4Alpha);
cCS = local_flatten_smoothing(context.C_CS,l4Alpha);
cCC = local_flatten_smoothing(context.C_CC,l4Alpha);
cCI = local_flatten_smoothing(context.C_CI,l4Alpha);
cIS = local_flatten_smoothing(context.C_IS,l4Alpha);
cIC = local_flatten_smoothing(context.C_IC,l4Alpha);
cII = local_flatten_smoothing(context.C_II,l4Alpha);
kControl = (1-l6Alpha)*kBase+l6Alpha*speye(n);

dL4ESS = (cSS*dSuseS+cSC*dCuseS)/context.L4SEp;
dL4ESC = (cSS*dSuseC+cSC*dCuseC)/context.L4SEp;
dL4ECS = (cCS*dSuseS+cCC*dCuseS)/context.L4CEp;
dL4ECC = (cCS*dSuseC+cCC*dCuseC)/context.L4CEp;
dL4EIS = (cIS*dSuseS+cIC*dCuseS)/context.L4IEp;
dL4EIC = (cIS*dSuseC+cIC*dCuseC)/context.L4IEp;
dL4ISI = (cSI*dIuseI)/context.L4SIp;
dL4ICI = (cCI*dIuseI)/context.L4CIp;
dL4III = (cII*dIuseI)/context.L4IIp;

dEuseS = wS*dSuseS+wC*dCuseS;
dEuseC = wS*dSuseC+wC*dCuseC;
dL6S = dL6Base*kControl*dEuseS;
dL6C = dL6Base*kControl*dEuseC;

DS4E = spdiags(dS4E,0,n,n); DS4I = spdiags(dS4I,0,n,n); DS6 = spdiags(dS6,0,n,n);
DC4E = spdiags(dC4E,0,n,n); DC4I = spdiags(dC4I,0,n,n); DC6 = spdiags(dC6,0,n,n);
DI4E = spdiags(dI4E,0,n,n); DI4I = spdiags(dI4I,0,n,n); DI6 = spdiags(dI6,0,n,n);

J = [DS4E*dL4ESS+DS6*dL6S, DS4E*dL4ESC+DS6*dL6C, DS4I*dL4ISI; ...
     DC4E*dL4ECS+DC6*dL6S, DC4E*dL4ECC+DC6*dL6C, DC4I*dL4ICI; ...
     DI4E*dL4EIS+DI6*dL6S, DI4E*dL4EIC+DI6*dL6C, DI4I*dL4III];
end

function controlled = local_flatten_smoothing(matrix,alpha)
n = size(matrix,1);
rowSums = full(sum(matrix,2));
if max(rowSums)-min(rowSums) > 1e-10*max(1,max(abs(rowSums)))
    error('SpatialControl:RowSums','Expected a translation-invariant block.');
end
controlled = (1-alpha)*matrix+alpha*mean(rowSums)*speye(n);
end

function [sUse,cUse,iUse,dSuseS,dSuseC,dCuseS,dCuseC,dIuseI] = ...
        local_preprocess_with_derivative(s,c,i,context)
n = numel(s);
if ~context.Isaturation
    sUse=s; cUse=c; iUse=i;
    identity=speye(n);
    dSuseS=identity; dSuseC=sparse(n,n);
    dCuseS=sparse(n,n); dCuseC=identity; dIuseI=identity;
    return
end
wS = 1-context.CWeight;
wC = context.CWeight;
eRaw = wS*s+wC*c;
eBase = L6Convert(eRaw,context.EKpUse);
adjustment = eBase./eRaw;
adjustment(~isfinite(adjustment))=1;
dEbase = L6Convert_grad(eRaw,context.EKpUse);
dAdjustment = (dEbase.*eRaw-eBase)./(eRaw.^2);
dAdjustment(~isfinite(dAdjustment))=0;
sUse=s.*adjustment;
cUse=c.*adjustment;
iUse=InhMulp(i,context.IKpUse);
dI = local_inhibition_gradient(i,context.IKpUse);
dSuseS=spdiags(adjustment+s.*dAdjustment*wS,0,n,n);
dSuseC=spdiags(s.*dAdjustment*wC,0,n,n);
dCuseS=spdiags(c.*dAdjustment*wS,0,n,n);
dCuseC=spdiags(adjustment+c.*dAdjustment*wC,0,n,n);
dIuseI=spdiags(dI,0,n,n);
end

function derivative = local_inhibition_gradient(x,parameters)
h = 1e-5*max(1,abs(x));
derivative = (InhMulp(x+h,parameters)-InhMulp(x-h,parameters))./(2*h);
derivative(~isfinite(derivative))=0;
end
