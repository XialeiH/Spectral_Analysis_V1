function J = l6ns_state_jacobian(setup,pathwayName,state,w)
% Analytic D Phi at an arbitrary equilibrium for the FPP interventions.

context = setup.Context;
state = state(:);
n = numel(state)/3;
if n~=prod(context.MapSize)
    error('State size is inconsistent with context.MapSize.');
end
s = state(1:n);
c = state(n+(1:n));
i = state(2*n+(1:n));

if strcmpi(pathwayName,'L6')
    l6FreezeWeight = w;
    inhibitoryGain = 1;
elseif any(strcmpi(pathwayName,{'I','Inhibition'}))
    l6FreezeWeight = 0;
    inhibitoryGain = 1-w;
else
    error('pathwayName must be L6, I, or Inhibition.');
end

[sUse,cUse,iUse,dSuseS,dSuseC,dCuseS,dCuseC,dIuseI] = ...
    local_preprocess_with_derivative(s,c,i,context);

l4ES = (context.C_SS*sUse + context.C_SC*cUse)/context.L4SEp;
l4EC = (context.C_CS*sUse + context.C_CC*cUse)/context.L4CEp;
l4EI = (context.C_IS*sUse + context.C_IC*cUse)/context.L4IEp;
l4IS = (context.C_SI*iUse)/context.L4SIp;
l4IC = (context.C_CI*iUse)/context.L4CIp;
l4II = (context.C_II*iUse)/context.L4IIp;

dL4ESS = (context.C_SS*dSuseS + context.C_SC*dCuseS)/context.L4SEp;
dL4ESC = (context.C_SS*dSuseC + context.C_SC*dCuseC)/context.L4SEp;
dL4ECS = (context.C_CS*dSuseS + context.C_CC*dCuseS)/context.L4CEp;
dL4ECC = (context.C_CS*dSuseC + context.C_CC*dCuseC)/context.L4CEp;
dL4EIS = (context.C_IS*dSuseS + context.C_IC*dCuseS)/context.L4IEp;
dL4EIC = (context.C_IS*dSuseC + context.C_IC*dCuseC)/context.L4IEp;
dL4ISI = inhibitoryGain*(context.C_SI*dIuseI)/context.L4SIp;
dL4ICI = inhibitoryGain*(context.C_CI*dIuseI)/context.L4CIp;
dL4III = inhibitoryGain*(context.C_II*dIuseI)/context.L4IIp;

if inhibitoryGain~=1
    fixed = context.FixedPoint(:);
    [~,~,iFixedUse] = local_preprocess_with_derivative( ...
        fixed(1:n),fixed(n+(1:n)),fixed(2*n+(1:n)),context);
    l4ISFixed = (context.C_SI*iFixedUse)/context.L4SIp;
    l4ICFixed = (context.C_CI*iFixedUse)/context.L4CIp;
    l4IIFixed = (context.C_II*iFixedUse)/context.L4IIp;
    l4IS = l4ISFixed + inhibitoryGain*(l4IS-l4ISFixed);
    l4IC = l4ICFixed + inhibitoryGain*(l4IC-l4ICFixed);
    l4II = l4IIFixed + inhibitoryGain*(l4II-l4IIFixed);
end

wS = 1-context.CWeight;
wC = context.CWeight;
eUse = wS*sUse+wC*cUse;
K = build_conv_matrix_circular(context.L6Kernel,context.MapSize(1),context.MapSize(2));
convolved = K*eUse;
l6Raw = L6Convert(reshape(convolved,context.MapSize),context.L6Parameters);
l6Raw = l6Raw(:)/3;
l6Dynamic = min(max(l6Raw,1),40);
l6Used = l6FreezeWeight*context.FixedL6LibInd(:) + ...
    (1-l6FreezeWeight)*l6Dynamic;

l6Gradient = L6Convert_grad(reshape(convolved,context.MapSize),context.L6Parameters);
l6Gradient = l6Gradient(:);
l6Mask = l6Raw>1 & l6Raw<40;
dL6Base = spdiags((l6Mask/3).*l6Gradient,0,n,n);
dEuseS = wS*dSuseS+wC*dCuseS;
dEuseC = wS*dSuseC+wC*dCuseC;
dL6S = (1-l6FreezeWeight)*(dL6Base*K*dEuseS);
dL6C = (1-l6FreezeWeight)*(dL6Base*K*dEuseC);

[dS4E,dS4I,dS6] = local_h96baseline_pref6D_grads('S',l4ES,l4IS,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
[dC4E,dC4I,dC6] = local_h96baseline_pref6D_grads('C',l4EC,l4IC,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
[dI4E,dI4I,dI6] = local_h96baseline_pref6D_grads('I',l4EI,l4II,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);

DS4E = spdiags(dS4E,0,n,n); DS4I = spdiags(dS4I,0,n,n); DS6 = spdiags(dS6,0,n,n);
DC4E = spdiags(dC4E,0,n,n); DC4I = spdiags(dC4I,0,n,n); DC6 = spdiags(dC6,0,n,n);
DI4E = spdiags(dI4E,0,n,n); DI4I = spdiags(dI4I,0,n,n); DI6 = spdiags(dI6,0,n,n);

J = [DS4E*dL4ESS+DS6*dL6S, DS4E*dL4ESC+DS6*dL6C, DS4I*dL4ISI; ...
     DC4E*dL4ECS+DC6*dL6S, DC4E*dL4ECC+DC6*dL6C, DC4I*dL4ICI; ...
     DI4E*dL4EIS+DI6*dL6S, DI4E*dL4EIC+DI6*dL6C, DI4I*dL4III];
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
