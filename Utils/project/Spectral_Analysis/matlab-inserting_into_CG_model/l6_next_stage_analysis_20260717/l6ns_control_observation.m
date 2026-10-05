function operators = l6ns_control_observation(context,cfg)
% Build proposal-specific stimulus controls and local/laminar observations.

fixed = context.FixedPoint(:);
n = numel(fixed)/3;
base = l6ns_phi(fixed,0,context,[1 1],1);
tau = cfg.TauMs;

contrastStep = cfg.ContrastProbeStep;
plusContext = context;
minusContext = context;
plusContext.ContrastUse = context.ContrastUse + contrastStep;
minusContext.ContrastUse = max(context.ContrastUse-contrastStep,eps);
bContrastPlus = (l6ns_phi(fixed,0,plusContext,[1 1],1)-base) / ...
    (contrastStep*tau);
bContrastMinus = (l6ns_phi(fixed,0,minusContext,[1 1],1)-base) / ...
    (contrastStep*tau);

orientationStep = cfg.OrientationProbeStepDeg;
plusContext = context;
minusContext = context;
plusContext.OrientationUse = context.OrientationUse + orientationStep;
minusContext.OrientationUse = context.OrientationUse - orientationStep;
bOrientationPlus = (l6ns_phi(fixed,0,plusContext,[1 1],1)-base) / ...
    (orientationStep*tau);
bOrientationMinus = (l6ns_phi(fixed,0,minusContext,[1 1],1)-base) / ...
    (orientationStep*tau);

mapSize = context.MapSize;
[yy,xx] = ndgrid(linspace(-1,1,mapSize(1)),linspace(-1,1,mapSize(2)));
surroundMask = sqrt(xx.^2+yy.^2) >= 0.55;
surroundStep = cfg.SurroundProbeStep;
contrastBase = context.ContrastUse;
if isscalar(contrastBase)
    contrastBase = repmat(contrastBase,n,1);
else
    contrastBase = contrastBase(:);
end
plusContext = context;
minusContext = context;
plusContext.ContrastUse = contrastBase + surroundStep*surroundMask(:);
minusContext.ContrastUse = max(contrastBase-surroundStep*surroundMask(:),eps);
bSurroundPlus = (l6ns_phi(fixed,0,plusContext,[1 1],1)-base) / ...
    (surroundStep*tau);
bSurroundMinus = (l6ns_phi(fixed,0,minusContext,[1 1],1)-base) / ...
    (surroundStep*tau);

b = [bContrastPlus,bContrastMinus,bOrientationPlus,bOrientationMinus, ...
    bSurroundPlus,bSurroundMinus];
controlNames = {'contrast+','contrast-','orientation+','orientation-', ...
    'surround+','surround-'};

s = fixed(1:n);
c = fixed(n+(1:n));
i = fixed(2*n+(1:n));
e = (1-context.CWeight)*s + context.CWeight*c;
[~,eIndex] = max(e);
[~,iIndex] = max(i);

cObs = sparse(5,3*n);
cObs(1,eIndex) = 1-context.CWeight;
cObs(1,n+eIndex) = context.CWeight;
cObs(2,2*n+iIndex) = 1;

[row,col] = ind2sub(mapSize,eIndex);
rowRange = max(1,row-1):min(mapSize(1),row+1);
colRange = max(1,col-1):min(mapSize(2),col+1);
localMask = false(mapSize);
localMask(rowRange,colRange) = true;
localIndex = find(localMask(:));
cObs(3,localIndex) = (1-context.CWeight)/numel(localIndex);
cObs(3,n+localIndex) = context.CWeight/numel(localIndex);
cObs(4,1:n) = (1-context.CWeight)/n;
cObs(4,n+(1:n)) = context.CWeight/n;
cObs(5,2*n+(1:n)) = 1/n;
observationNames = {'local_E_unit','local_I_unit','local_E_MUA', ...
    'laminar_E_mean','laminar_I_mean'};

operators = struct();
operators.B = sparse(b);
operators.C = cObs;
operators.ControlNames = controlNames;
operators.ObservationNames = observationNames;
operators.ControlNorms = vecnorm(b)';
operators.LocalEIndex = eIndex;
operators.LocalIIndex = iIndex;
operators.LocalEMapSubscript = [row col];
operators.SurroundMask = surroundMask;
end
