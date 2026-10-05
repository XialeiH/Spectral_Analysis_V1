function J=l6ns_controlled_jacobian(setup,operators,fixedPoint)
% Analytic D Phi for a controlled model at that model's own equilibrium.

context=setup.Context;
state=fixedPoint(:);
n=numel(state)/3;
s=state(1:n); c=state(n+(1:n)); i=state(2*n+(1:n));
[sUse,cUse,iUse,dSuseS,dSuseC,dCuseS,dCuseC,dIuseI]= ...
    local_preprocess_with_derivative(s,c,i,context);

l4ES=(operators.C_SS*sUse+operators.C_SC*cUse)/context.L4SEp;
l4EC=(operators.C_CS*sUse+operators.C_CC*cUse)/context.L4CEp;
l4EI=(operators.C_IS*sUse+operators.C_IC*cUse)/context.L4IEp;
l4IS=(operators.C_SI*iUse)/context.L4SIp;
l4IC=(operators.C_CI*iUse)/context.L4CIp;
l4II=(operators.C_II*iUse)/context.L4IIp;

wS=1-context.CWeight; wC=context.CWeight;
eUse=wS*sUse+wC*cUse;
convolved=operators.L6*eUse;
l6Raw=L6Convert(reshape(convolved,context.MapSize),context.L6Parameters);
l6Raw=l6Raw(:)/3;
l6Used=min(max(l6Raw,1),40);
l6Gradient=L6Convert_grad(reshape(convolved,context.MapSize), ...
    context.L6Parameters);
l6Mask=l6Raw>1 & l6Raw<40;
dL6Base=spdiags((l6Mask/3).*l6Gradient(:),0,n,n);

[dS4E,dS4I,dS6]=local_h96baseline_pref6D_grads('S',l4ES,l4IS,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
[dC4E,dC4I,dC6]=local_h96baseline_pref6D_grads('C',l4EC,l4IC,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
[dI4E,dI4I,dI6]=local_h96baseline_pref6D_grads('I',l4EI,l4II,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);

cSS=operators.C_SS; cSC=operators.C_SC; cSI=operators.C_SI;
cCS=operators.C_CS; cCC=operators.C_CC; cCI=operators.C_CI;
cIS=operators.C_IS; cIC=operators.C_IC; cII=operators.C_II;
dL4ESS=(cSS*dSuseS+cSC*dCuseS)/context.L4SEp;
dL4ESC=(cSS*dSuseC+cSC*dCuseC)/context.L4SEp;
dL4ECS=(cCS*dSuseS+cCC*dCuseS)/context.L4CEp;
dL4ECC=(cCS*dSuseC+cCC*dCuseC)/context.L4CEp;
dL4EIS=(cIS*dSuseS+cIC*dCuseS)/context.L4IEp;
dL4EIC=(cIS*dSuseC+cIC*dCuseC)/context.L4IEp;
dL4ISI=(cSI*dIuseI)/context.L4SIp;
dL4ICI=(cCI*dIuseI)/context.L4CIp;
dL4III=(cII*dIuseI)/context.L4IIp;
dEuseS=wS*dSuseS+wC*dCuseS;
dEuseC=wS*dSuseC+wC*dCuseC;
dL6S=dL6Base*operators.L6*dEuseS;
dL6C=dL6Base*operators.L6*dEuseC;

DS4E=spdiags(dS4E,0,n,n); DS4I=spdiags(dS4I,0,n,n); DS6=spdiags(dS6,0,n,n);
DC4E=spdiags(dC4E,0,n,n); DC4I=spdiags(dC4I,0,n,n); DC6=spdiags(dC6,0,n,n);
DI4E=spdiags(dI4E,0,n,n); DI4I=spdiags(dI4I,0,n,n); DI6=spdiags(dI6,0,n,n);
J=[DS4E*dL4ESS+DS6*dL6S, DS4E*dL4ESC+DS6*dL6C, DS4I*dL4ISI; ...
   DC4E*dL4ECS+DC6*dL6S, DC4E*dL4ECC+DC6*dL6C, DC4I*dL4ICI; ...
   DI4E*dL4EIS+DI6*dL6S, DI4E*dL4EIC+DI6*dL6C, DI4I*dL4III];
end

function [sUse,cUse,iUse,dSuseS,dSuseC,dCuseS,dCuseC,dIuseI]= ...
        local_preprocess_with_derivative(s,c,i,context)
n=numel(s);
if ~context.Isaturation
    sUse=s; cUse=c; iUse=i; identity=speye(n);
    dSuseS=identity; dSuseC=sparse(n,n); dCuseS=sparse(n,n);
    dCuseC=identity; dIuseI=identity; return
end
wS=1-context.CWeight; wC=context.CWeight;
eRaw=wS*s+wC*c;
eBase=L6Convert(eRaw,context.EKpUse);
adjustment=eBase./eRaw; adjustment(~isfinite(adjustment))=1;
dEbase=L6Convert_grad(eRaw,context.EKpUse);
dAdjustment=(dEbase.*eRaw-eBase)./(eRaw.^2);
dAdjustment(~isfinite(dAdjustment))=0;
sUse=s.*adjustment; cUse=c.*adjustment; iUse=InhMulp(i,context.IKpUse);
h=1e-5*max(1,abs(i));
dI=(InhMulp(i+h,context.IKpUse)-InhMulp(i-h,context.IKpUse))./(2*h);
dI(~isfinite(dI))=0;
dSuseS=spdiags(adjustment+s.*dAdjustment*wS,0,n,n);
dSuseC=spdiags(s.*dAdjustment*wC,0,n,n);
dCuseS=spdiags(c.*dAdjustment*wS,0,n,n);
dCuseC=spdiags(adjustment+c.*dAdjustment*wC,0,n,n);
dIuseI=spdiags(dI,0,n,n);
end
