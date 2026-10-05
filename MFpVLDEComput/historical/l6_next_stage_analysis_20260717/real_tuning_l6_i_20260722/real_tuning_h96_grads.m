function [dPhi_dL4E,dPhi_dL4I,dPhi_dL6] = real_tuning_h96_grads( ...
    celltype,L4EUse,L4IUse,L6Use,ContrastUse,OrientationUse,PixLGNCtgr)
% Derivatives of the exact domain-extended h96 response used by Phi.

LGNnum = size(PixLGNCtgr,2);
dPhi_dL4E = zeros(size(L4EUse));
dPhi_dL4I = zeros(size(L4IUse));
dPhi_dL6 = zeros(size(L6Use));
thetaMap = [0 135 90 45];

for LGNInd = 1:LGNnum
    contrastEff = ContrastUse;
    oriCos2 = ones(size(OrientationUse));
    if LGNInd <= 4
        gammaUse = abs(local_wrap_to_90(double(OrientationUse)-thetaMap(LGNInd)));
        oriCos2 = cosd(2*gammaUse);
    else
        contrastEff(:) = 0;
    end
    [~,gE,gI,gL6] = h96_pref6D_domain_extended_response( ...
        celltype,L4EUse,L4IUse,L6Use,contrastEff,oriCos2);
    outsideTrainingBox = L4EUse>47500 | L4IUse>35315.63671875;
    if any(outsideTrainingBox,'all')
        [gEExtended,gIExtended,gL6Extended] = local_extended_derivatives( ...
            celltype,L4EUse,L4IUse,L6Use,contrastEff,oriCos2);
        gE(outsideTrainingBox) = gEExtended(outsideTrainingBox);
        gI(outsideTrainingBox) = gIExtended(outsideTrainingBox);
        gL6(outsideTrainingBox) = gL6Extended(outsideTrainingBox);
    end
    weights = PixLGNCtgr(:,LGNInd);
    dPhi_dL4E = dPhi_dL4E + weights.*gE;
    dPhi_dL4I = dPhi_dL4I + weights.*gI;
    dPhi_dL6 = dPhi_dL6 + weights.*gL6;
end
end

function [gE,gI,gL6] = local_extended_derivatives( ...
    celltype,L4EUse,L4IUse,L6Use,contrastEff,oriCos2)
% The extension slope depends on boundary-network gradients, so its exact
% first derivative would require NN Hessians. Differentiate that local,
% elementwise extension numerically while retaining analytic NN gradients
% everywhere inside the trained box.
hE = 1e-5*max(1,abs(double(L4EUse)));
hI = 1e-5*max(1,abs(double(L4IUse)));
h6 = 1e-5*max(1,abs(double(L6Use)));
gE = (h96_pref6D_domain_extended_response(celltype,L4EUse+hE,L4IUse, ...
    L6Use,contrastEff,oriCos2)- ...
    h96_pref6D_domain_extended_response(celltype,L4EUse-hE,L4IUse, ...
    L6Use,contrastEff,oriCos2))./(2*hE);
gI = (h96_pref6D_domain_extended_response(celltype,L4EUse,L4IUse+hI, ...
    L6Use,contrastEff,oriCos2)- ...
    h96_pref6D_domain_extended_response(celltype,L4EUse,L4IUse-hI, ...
    L6Use,contrastEff,oriCos2))./(2*hI);
gL6 = (h96_pref6D_domain_extended_response(celltype,L4EUse,L4IUse, ...
    L6Use+h6,contrastEff,oriCos2)- ...
    h96_pref6D_domain_extended_response(celltype,L4EUse,L4IUse, ...
    L6Use-h6,contrastEff,oriCos2))./(2*h6);
end

function thetaWrapped = local_wrap_to_90(thetaDeg)
thetaWrapped = mod(thetaDeg+90,180)-90;
end
