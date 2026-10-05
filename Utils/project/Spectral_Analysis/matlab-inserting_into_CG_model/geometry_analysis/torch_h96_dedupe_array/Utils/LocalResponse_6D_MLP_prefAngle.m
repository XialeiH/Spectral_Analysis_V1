function LDEOutLIBy = LocalResponse_6D_MLP_prefAngle(celltype, alphaUse, LGNctgrUse, L4EUse, L4IUse, L6Use, ContrastUse) %#codegen
% Translate runtime variables into the 5D cos-only training contract:
% [L4E, L4I, L6, ContrastEff, OriCos2].
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
        LDEOutLIBy = predict_pref6D_S(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
    case 'C'
        LDEOutLIBy = predict_pref6D_C(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
    case 'I'
        LDEOutLIBy = predict_pref6D_I(L4EUse, L4IUse, L6Use, contrastEff, oriCos2);
    otherwise
        error('Unknown celltype');
end
end

function thetaWrapped = wrapTo90_pref(thetaDeg)
thetaWrapped = mod(thetaDeg + 90, 180) - 90;
end
