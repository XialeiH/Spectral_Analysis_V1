function [LDEOutLIBy, dL4E, dL4I, dL6] = LocalResponse_6D_MLP_prefAngle_with_grad(celltype, alphaUse, LGNctgrUse, L4EUse, L4IUse, L6Use, ContrastUse) %#codegen
% State derivatives for the h96 baseline pref-angle NN.
if ~isequal(size(alphaUse), size(LGNctgrUse), size(L4EUse), size(L4IUse), size(L6Use), size(ContrastUse))
    error('All LocalResponse inputs must share shape.');
end

thetaPref = zeros(size(alphaUse));
gammaUse = zeros(size(alphaUse));
oriCos2 = ones(size(alphaUse));
contrastEff = ContrastUse;

fgMask = (LGNctgrUse >= 1) & (LGNctgrUse <= 4);
if any(fgMask(:))
    thetaMap = [0, 135, 90, 45];
    thetaPref(fgMask) = thetaMap(double(LGNctgrUse(fgMask)));
    gammaUse(fgMask) = abs(wrapTo90_pref(double(alphaUse(fgMask)) - thetaPref(fgMask)));
    oriCos2(fgMask) = cosd(2 .* gammaUse(fgMask));
end

bgMask = (LGNctgrUse == 5);
if any(bgMask(:))
    contrastEff(bgMask) = 0;
    gammaUse(bgMask) = abs(wrapTo90_pref(double(alphaUse(bgMask)) - thetaPref(bgMask)));
    oriCos2(bgMask) = cosd(2 .* gammaUse(bgMask));
end

switch celltype
    case 'S'
        [yParent, dL4E, dL4I, dL6] = predict_pref6D_S_with_grad(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
    case 'C'
        [yParent, dL4E, dL4I, dL6] = predict_pref6D_C_with_grad(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
    case 'I'
        [yParent, dL4E, dL4I, dL6] = predict_pref6D_I_with_grad(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
    otherwise
        error('Unknown celltype');
end

LDEOutLIBy = yParent;
end

function thetaWrapped = wrapTo90_pref(thetaDeg)
thetaWrapped = mod(thetaDeg + 90, 180) - 90;
end
