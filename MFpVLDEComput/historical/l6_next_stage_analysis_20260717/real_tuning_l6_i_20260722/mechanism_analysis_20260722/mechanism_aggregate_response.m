function [output,dL4E,dL4I,dL6] = mechanism_aggregate_response( ...
    celltype,L4E,L4I,L6,context,responseMode)
% Evaluate the weighted five-LGN response with or without domain extension.

if nargin<6
    responseMode = 'extended';
end
n = numel(L4E);
output = zeros(n,1);
dL4E = zeros(n,1);
dL4I = zeros(n,1);
dL6 = zeros(n,1);
thetaMap = [0 135 90 45];

for lgnIndex = 1:5
    contrastEff = context.ContrastUse;
    oriCos2 = ones(size(context.OrientationUse));
    if lgnIndex<=4
        gamma = abs(local_wrap_to_90(double(context.OrientationUse)-thetaMap(lgnIndex)));
        oriCos2 = cosd(2*gamma);
    else
        contrastEff(:) = 0;
    end
    if strcmp(responseMode,'extended')
        [y,gE,gI,g6] = h96_pref6D_domain_extended_response( ...
            celltype,L4E,L4I,L6,contrastEff,oriCos2);
    elseif strcmp(responseMode,'raw')
        [y,gE,gI,g6] = local_raw_predictor( ...
            celltype,L4E,L4I,L6,contrastEff,oriCos2);
    else
        error('Mechanism:ResponseMode','Unknown response mode %s.',responseMode);
    end
    weight = context.PixLGNCtgr(:,lgnIndex);
    output = output+weight.*y;
    dL4E = dL4E+weight.*gE;
    dL4I = dL4I+weight.*gI;
    dL6 = dL6+weight.*g6;
end
end

function [y,gE,gI,g6] = local_raw_predictor(celltype,L4E,L4I,L6,contrastEff,oriCos2)
switch celltype
    case 'S'
        [y,gE,gI,g6] = predict_pref6D_S_with_grad( ...
            L4E,L4I,L6,contrastEff,oriCos2);
    case 'C'
        [y,gE,gI,g6] = predict_pref6D_C_with_grad( ...
            L4E,L4I,L6,contrastEff,oriCos2);
    case 'I'
        [y,gE,gI,g6] = predict_pref6D_I_with_grad( ...
            L4E,L4I,L6,contrastEff,oriCos2);
    otherwise
        error('Mechanism:CellType','Unknown cell type %s.',celltype);
end
end

function wrapped = local_wrap_to_90(angleDeg)
wrapped = mod(angleDeg+90,180)-90;
end
