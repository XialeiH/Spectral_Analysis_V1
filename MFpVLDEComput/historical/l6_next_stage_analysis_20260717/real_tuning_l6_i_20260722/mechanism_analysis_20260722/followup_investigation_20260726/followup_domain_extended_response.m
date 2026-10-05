function y = followup_domain_extended_response(celltype,L4E,L4I,L6, ...
        contrastEff,oriCos2,edgeScale,blendScale)
% Evaluate the existing high-input extension with controlled edges and width.

edgeE = 47500*edgeScale;
edgeI = 35315.63671875*edgeScale;
blendE = 5000*blendScale;
blendI = 5000*blendScale;
rateHigh = 200;
minimumK = 1e-4;
yNN = local_predict(celltype,L4E,L4I,L6,contrastEff,oriCos2,false);
upper = L4E>edgeE | L4I>edgeI;
if ~any(upper(:)); y=yNN; return; end

boundaryE = min(L4E,edgeE);
boundaryI = min(L4I,edgeI);
yBoundary = local_predict(celltype,boundaryE,boundaryI,L6, ...
    contrastEff,oriCos2,false);
[dE,dI] = local_boundary_gradient(celltype,boundaryE,boundaryI,L6, ...
    contrastEff,oriCos2);
exE = max(0,double(L4E)-edgeE);
exI = max(0,double(L4I)-edgeI);
distance = sqrt((exE/blendE).^2+(exI/blendI).^2);
weight = ones(size(distance));
insideBlend = distance>0 & distance<1;
t = distance(insideBlend);
weight(insideBlend)=6*t.^5-15*t.^4+10*t.^3;
weight(distance<=0)=0;

radius = sqrt(exE.^2+exI.^2);
uE = zeros(size(radius)); uI = zeros(size(radius));
positive = radius>0;
uE(positive)=exE(positive)./radius(positive);
uI(positive)=exI(positive)./radius(positive);
outSlope=dE.*uE+dI.*uI;
ySig=yBoundary;
gap=max(rateHigh-yBoundary(upper),1e-6);
k=max(outSlope(upper)./gap,minimumK);
q=min(k.*radius(upper),50);
ySig(upper)=rateHigh-gap.*exp(-q);
y=yNN+weight.*(ySig-yNN);
end

function [dE,dI] = local_boundary_gradient(celltype,E,I,L6,contrast,orientation)
hE=1e-4*max(1,abs(double(E)));
hI=1e-4*max(1,abs(double(I)));
dE=(local_predict(celltype,E+hE,I,L6,contrast,orientation,false)- ...
    local_predict(celltype,E-hE,I,L6,contrast,orientation,false))./(2*hE);
dI=(local_predict(celltype,E,I+hI,L6,contrast,orientation,false)- ...
    local_predict(celltype,E,I-hI,L6,contrast,orientation,false))./(2*hI);
dE(~isfinite(dE))=0; dI(~isfinite(dI))=0;
end

function y = local_predict(celltype,E,I,L6,contrast,orientation,needGradient) %#ok<INUSD>
switch celltype
    case 'S'; y=predict_pref6D_S(E,I,L6,contrast,orientation);
    case 'C'; y=predict_pref6D_C(E,I,L6,contrast,orientation);
    case 'I'; y=predict_pref6D_I(E,I,L6,contrast,orientation);
    otherwise; error('Followup:CellType','Unknown cell type %s.',celltype);
end
end
