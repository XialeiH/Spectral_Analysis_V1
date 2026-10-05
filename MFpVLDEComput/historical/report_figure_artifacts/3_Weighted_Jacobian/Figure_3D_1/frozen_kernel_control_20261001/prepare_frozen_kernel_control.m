function validation = prepare_frozen_kernel_control
% Freeze all non-kernel derivatives at the original Figure 3D baseline state.
root = fileparts(mfilename('fullpath'));
runtime = fullfile(root,'runtime');
addpath(runtime,'-begin');
setupFile = fullfile(root,'input','global_bifurcation_setup.mat');
sourceFile = fullfile(fileparts(root),'source','baseline.mat');
loaded = load(setupFile,'setup');
context = loaded.setup.Context;
source = load(sourceFile,'fixedPoint','operators','J','eigenvalues');
state = source.fixedPoint(:);
baseline = source.operators;
n = numel(state)/3;
s = state(1:n); c = state(n+(1:n)); i = state(2*n+(1:n));
wS = 1-context.CWeight; wC = context.CWeight;
if context.Isaturation
    eRaw = wS*s+wC*c;
    eBase = L6Convert(eRaw,context.EKpUse);
    adjustment = eBase./eRaw; adjustment(~isfinite(adjustment)) = 1;
    dEbase = L6Convert_grad(eRaw,context.EKpUse);
    dAdjustment = (dEbase.*eRaw-eBase)./(eRaw.^2);
    dAdjustment(~isfinite(dAdjustment)) = 0;
    sUse = s.*adjustment; cUse = c.*adjustment;
    iUse = InhMulp(i,context.IKpUse);
    h = 1e-5*max(1,abs(i));
    dI = (InhMulp(i+h,context.IKpUse)-InhMulp(i-h,context.IKpUse))./(2*h);
    dI(~isfinite(dI)) = 0;
    d.SS = spdiags(adjustment+s.*dAdjustment*wS,0,n,n);
    d.SC = spdiags(s.*dAdjustment*wC,0,n,n);
    d.CS = spdiags(c.*dAdjustment*wS,0,n,n);
    d.CC = spdiags(adjustment+c.*dAdjustment*wC,0,n,n);
    d.II = spdiags(dI,0,n,n);
else
    sUse = s; cUse = c; iUse = i;
    d.SS = speye(n); d.SC = sparse(n,n); d.CS = sparse(n,n);
    d.CC = speye(n); d.II = speye(n);
end

l4ES = (baseline.C_SS*sUse+baseline.C_SC*cUse)/context.L4SEp;
l4EC = (baseline.C_CS*sUse+baseline.C_CC*cUse)/context.L4CEp;
l4EI = (baseline.C_IS*sUse+baseline.C_IC*cUse)/context.L4IEp;
l4IS = (baseline.C_SI*iUse)/context.L4SIp;
l4IC = (baseline.C_CI*iUse)/context.L4CIp;
l4II = (baseline.C_II*iUse)/context.L4IIp;
convolved = baseline.L6*(wS*sUse+wC*cUse);
l6Raw = L6Convert(reshape(convolved,context.MapSize),context.L6Parameters);
l6Raw = l6Raw(:)/3;
l6Used = min(max(l6Raw,1),40);
l6Gradient = L6Convert_grad(reshape(convolved,context.MapSize),context.L6Parameters);
l6Mask = l6Raw>1 & l6Raw<40;
[dS4E,dS4I,dS6] = local_h96baseline_pref6D_grads('S',l4ES,l4IS,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
[dC4E,dC4I,dC6] = local_h96baseline_pref6D_grads('C',l4EC,l4IC,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
[dI4E,dI4I,dI6] = local_h96baseline_pref6D_grads('I',l4EI,l4II,l6Used, ...
    context.ContrastUse,context.OrientationUse,context.PixLGNCtgr);
gradientNames = {'S4E','S4I','S6','C4E','C4I','C6','I4E','I4I','I6'};
gradientValues = {dS4E,dS4I,dS6,dC4E,dC4I,dC6,dI4E,dI4I,dI6};
for index = 1:numel(gradientNames)
    assert(all(isfinite(gradientValues{index})));
    response.(gradientNames{index}) = spdiags(gradientValues{index},0,n,n);
end
normalizationNames = {'L4SEp','L4CEp','L4IEp','L4SIp','L4CIp','L4IIp'};
for index = 1:numel(normalizationNames)
    normalization.(normalizationNames{index}) = context.(normalizationNames{index});
end
frozen = struct('FixedPoint',state,'Preprocess',d,'Response',response, ...
    'L6Derivative',spdiags((l6Mask/3).*l6Gradient(:),0,n,n), ...
    'L6Mask',l6Mask,'Normalization',normalization,'CWeight',context.CWeight, ...
    'MapSize',context.MapSize,'Isaturation',context.Isaturation, ...
    'Inputs',struct('L4ES',l4ES,'L4EC',l4EC,'L4EI',l4EI, ...
    'L4IS',l4IS,'L4IC',l4IC,'L4II',l4II,'L6Raw',l6Raw,'L6Used',l6Used), ...
    'SourceFile',sourceFile,'SetupFile',setupFile);
[J0,~,rowErrors] = frozen_kernel_jacobian(frozen,baseline,0,0);
validation = struct('RelativeFrobeniusError',norm(J0-source.J,'fro')/norm(source.J,'fro'), ...
    'MaximumEntryError',max(abs(nonzeros(J0-source.J))), ...
    'MaximumRowSumError',max(rowErrors),'N',size(J0,1));
if isempty(validation.MaximumEntryError); validation.MaximumEntryError = 0; end
assert(validation.RelativeFrobeniusError < 1e-10, ...
    'Reconstructed baseline must match the original Figure 3D Jacobian.');
assert(validation.MaximumEntryError < 1e-8);
baselineEigenvalues = source.eigenvalues;
save(fullfile(root,'frozen_derivatives.mat'),'frozen','baseline','J0', ...
    'baselineEigenvalues','validation','-v7.3');
disp(validation);
fprintf('Saved derivatives once at the original baseline; reconstruction gate passed.\n');
end
