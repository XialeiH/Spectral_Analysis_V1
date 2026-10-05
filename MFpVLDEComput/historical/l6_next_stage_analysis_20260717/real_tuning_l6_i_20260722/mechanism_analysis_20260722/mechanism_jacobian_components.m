function [J,components,details] = mechanism_jacobian_components( ...
    state,context,gain6,gainI,responseMode)
% Exact true-tuning Jacobian split by recurrent input pathway.

if nargin<5
    responseMode = 'extended';
end

state = state(:);
n = numel(state)/3;
s = state(1:n);
c = state(n+(1:n));
i = state(2*n+(1:n));
wC = context.CWeight;
wS = 1-wC;

if context.Isaturation
    eRaw = wS*s+wC*c;
    eBase = L6Convert(eRaw,context.EKpUse);
    eAdjustment = eBase./eRaw;
    dEbase = L6Convert_grad(eRaw,context.EKpUse);
    dAdjustment = (dEbase.*eRaw-eBase)./(eRaw.^2);
    eAdjustment(~isfinite(eAdjustment)) = 1;
    dAdjustment(~isfinite(dAdjustment)) = 0;
    sUse = s.*eAdjustment;
    cUse = c.*eAdjustment;
    iUse = InhMulp(i,context.IKpUse);
    dIuse = local_inhibition_gradient(i,context.IKpUse);
else
    eAdjustment = ones(n,1);
    dAdjustment = zeros(n,1);
    sUse = s;
    cUse = c;
    iUse = i;
    dIuse = ones(n,1);
end

dSuseS = spdiags(eAdjustment+s.*dAdjustment*wS,0,n,n);
dSuseC = spdiags(s.*dAdjustment*wC,0,n,n);
dCuseS = spdiags(c.*dAdjustment*wS,0,n,n);
dCuseC = spdiags(eAdjustment+c.*dAdjustment*wC,0,n,n);
dIuseI = spdiags(dIuse,0,n,n);

l4ES = (context.C_SS*sUse+context.C_SC*cUse)/context.L4SEp;
l4EC = (context.C_CS*sUse+context.C_CC*cUse)/context.L4CEp;
l4EI = (context.C_IS*sUse+context.C_IC*cUse)/context.L4IEp;
l4IS = gainI*(context.C_SI*iUse)/context.L4SIp;
l4IC = gainI*(context.C_CI*iUse)/context.L4CIp;
l4II = gainI*(context.C_II*iUse)/context.L4IIp;

dL4ESS = (context.C_SS*dSuseS+context.C_SC*dCuseS)/context.L4SEp;
dL4ESC = (context.C_SS*dSuseC+context.C_SC*dCuseC)/context.L4SEp;
dL4ECS = (context.C_CS*dSuseS+context.C_CC*dCuseS)/context.L4CEp;
dL4ECC = (context.C_CS*dSuseC+context.C_CC*dCuseC)/context.L4CEp;
dL4EIS = (context.C_IS*dSuseS+context.C_IC*dCuseS)/context.L4IEp;
dL4EIC = (context.C_IS*dSuseC+context.C_IC*dCuseC)/context.L4IEp;
dL4ISI = gainI*(context.C_SI*dIuseI)/context.L4SIp;
dL4ICI = gainI*(context.C_CI*dIuseI)/context.L4CIp;
dL4III = gainI*(context.C_II*dIuseI)/context.L4IIp;

eUse = wS*sUse+wC*cUse;
K = build_conv_matrix_circular(context.L6Kernel,context.MapSize(1),context.MapSize(2));
convolved = K*eUse;
l6Unclamped = L6Convert(reshape(convolved,context.MapSize),context.L6Parameters);
l6Unclamped = l6Unclamped(:)/3;
l6Base = min(max(l6Unclamped,1),40);
l6Input = gain6*l6Base;
l6Mask = l6Unclamped>1 & l6Unclamped<40;
l6CurveGradient = L6Convert_grad(reshape(convolved,context.MapSize), ...
    context.L6Parameters);
l6CurveGradient = l6CurveGradient(:);
dL6base = spdiags((l6Mask/3).*l6CurveGradient,0,n,n);
dEuseS = wS*dSuseS+wC*dCuseS;
dEuseC = wS*dSuseC+wC*dCuseC;
dL6S = gain6*dL6base*K*dEuseS;
dL6C = gain6*dL6base*K*dEuseC;

if strcmp(responseMode,'extended')
    [sE,sI,s6] = real_tuning_h96_grads('S',l4ES,l4IS,l6Input, ...
        context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
    [cE,cI,c6] = real_tuning_h96_grads('C',l4EC,l4IC,l6Input, ...
        context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
    [iE,iI,i6] = real_tuning_h96_grads('I',l4EI,l4II,l6Input, ...
        context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
elseif strcmp(responseMode,'raw')
    [~,sE,sI,s6] = mechanism_aggregate_response( ...
        'S',l4ES,l4IS,l6Input,context,'raw');
    [~,cE,cI,c6] = mechanism_aggregate_response( ...
        'C',l4EC,l4IC,l6Input,context,'raw');
    [~,iE,iI,i6] = mechanism_aggregate_response( ...
        'I',l4EI,l4II,l6Input,context,'raw');
else
    error('Mechanism:ResponseMode','Unknown response mode %s.',responseMode);
end

sEdiag = spdiags(sE,0,n,n); sIdiag = spdiags(sI,0,n,n); s6diag = spdiags(s6,0,n,n);
cEdiag = spdiags(cE,0,n,n); cIdiag = spdiags(cI,0,n,n); c6diag = spdiags(c6,0,n,n);
iEdiag = spdiags(iE,0,n,n); iIdiag = spdiags(iI,0,n,n); i6diag = spdiags(i6,0,n,n);
zero = sparse(n,n);

components.Excitatory = [sEdiag*dL4ESS,sEdiag*dL4ESC,zero; ...
    cEdiag*dL4ECS,cEdiag*dL4ECC,zero; ...
    iEdiag*dL4EIS,iEdiag*dL4EIC,zero];
components.L6 = [s6diag*dL6S,s6diag*dL6C,zero; ...
    c6diag*dL6S,c6diag*dL6C,zero; ...
    i6diag*dL6S,i6diag*dL6C,zero];
components.Inhibition = [zero,zero,sIdiag*dL4ISI; ...
    zero,zero,cIdiag*dL4ICI; ...
    zero,zero,iIdiag*dL4III];
J = sparse(components.Excitatory+components.L6+components.Inhibition);

details.EAdjustment = eAdjustment;
details.IGradient = dIuse;
details.L4EGradient = [sE;cE;iE];
details.L4IGradient = [sI;cI;iI];
details.L6Gradient = [s6;c6;i6];
details.L4E = [l4ES;l4EC;l4EI];
details.L4I = [l4IS;l4IC;l4II];
details.L6 = l6Input;
end

function derivative = local_inhibition_gradient(values,parameters)
step = 1e-5*max(1,abs(values));
derivative = (InhMulp(values+step,parameters)- ...
    InhMulp(values-step,parameters))./(2*step);
derivative(~isfinite(derivative)) = 0;
end
