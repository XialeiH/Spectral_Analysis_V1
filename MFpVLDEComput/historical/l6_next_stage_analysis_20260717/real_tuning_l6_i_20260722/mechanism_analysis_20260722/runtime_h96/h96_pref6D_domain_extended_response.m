function [y, dL4E, dL4I, dL6] = h96_pref6D_domain_extended_response(celltype, L4EUse, L4IUse, L6Use, ContrastEff, OriCos2) %#codegen
% Extend the h96 pref6D response above the trained L4E/L4I box.
% Inside the trained box, this returns the original h96 predictor exactly.
L4EMax = 47500.0;
L4IMax = 35315.63671875;
blendE = 5000.0;
blendI = 5000.0;
rateHigh = 200.0;
minSigmoidRate = 1e-4;

upperOut = (L4EUse > L4EMax) | (L4IUse > L4IMax);
if ~any(upperOut(:))
    if nargout <= 1
        y = call_h96_pref6D_predictor(celltype, L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, false);
        dL4E = zeros(size(y));
        dL4I = zeros(size(y));
        dL6 = zeros(size(y));
    else
        [y, dL4E, dL4I, dL6] = call_h96_pref6D_predictor(celltype, L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, true);
    end
    return
end

needGrad = nargout > 1;
if needGrad
    [yNN, dE_NN, dI_NN, dL6_NN] = call_h96_pref6D_predictor(celltype, L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, true);
else
    yNN = call_h96_pref6D_predictor(celltype, L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, false);
    dE_NN = zeros(size(yNN));
    dI_NN = zeros(size(yNN));
    dL6_NN = zeros(size(yNN));
end

L4EBdry = min(L4EUse, L4EMax);
L4IBdry = min(L4IUse, L4IMax);
if needGrad
    [yBdry, dE_Bdry, dI_Bdry, dL6_Bdry] = call_h96_pref6D_predictor(celltype, L4EBdry, L4IBdry, L6Use, ContrastEff, OriCos2, true);
else
    yBdry = call_h96_pref6D_predictor(celltype, L4EBdry, L4IBdry, L6Use, ContrastEff, OriCos2, false);
    [dE_Bdry, dI_Bdry] = boundary_directional_finite_diff(celltype, L4EBdry, L4IBdry, L6Use, ContrastEff, OriCos2);
    dL6_Bdry = zeros(size(yBdry));
end

exE = max(0.0, double(L4EUse) - L4EMax);
exI = max(0.0, double(L4IUse) - L4IMax);
exEn = exE ./ blendE;
exIn = exI ./ blendI;
dBlend = sqrt(exEn.^2 + exIn.^2);

blendW = ones(size(dBlend));
blendDW = zeros(size(dBlend));
blendMask = dBlend > 0 & dBlend < 1;
t = dBlend(blendMask);
blendW(blendMask) = 6 .* t.^5 - 15 .* t.^4 + 10 .* t.^3;
blendDW(blendMask) = 30 .* t.^4 - 60 .* t.^3 + 30 .* t.^2;
blendW(dBlend <= 0) = 0;

dd_dE = zeros(size(dBlend));
dd_dI = zeros(size(dBlend));
posD = dBlend > 0;
dd_dE(posD) = exEn(posD) ./ (dBlend(posD) .* blendE);
dd_dI(posD) = exIn(posD) ./ (dBlend(posD) .* blendI);
dW_dE = blendDW .* dd_dE;
dW_dI = blendDW .* dd_dI;

rPhys = sqrt(exE.^2 + exI.^2);
uE = zeros(size(rPhys));
uI = zeros(size(rPhys));
posR = rPhys > 0;
uE(posR) = exE(posR) ./ rPhys(posR);
uI(posR) = exI(posR) ./ rPhys(posR);

outSlope = dE_Bdry .* uE + dI_Bdry .* uI;

ySig = yBdry;
dE_Sig = zeros(size(yBdry));
dI_Sig = zeros(size(yBdry));
dL6_Sig = dL6_Bdry;

highMask = upperOut;
if any(highMask(:))
    gap = max(rateHigh - yBdry(highMask), 1e-6);
    k = max(outSlope(highMask) ./ gap, minSigmoidRate);
    q = min(k .* rPhys(highMask), 50.0);
    eq = exp(-q);
    ySig(highMask) = rateHigh - gap .* eq;
    dyDr = gap .* k .* eq;
    dE_Sig(highMask) = dyDr .* uE(highMask);
    dI_Sig(highMask) = dyDr .* uI(highMask);
    dL6_Sig(highMask) = eq .* dL6_Bdry(highMask);
end

delta = ySig - yNN;
y = yNN + blendW .* delta;
dL4E = dE_NN + blendW .* (dE_Sig - dE_NN) + dW_dE .* delta;
dL4I = dI_NN + blendW .* (dI_Sig - dI_NN) + dW_dI .* delta;
dL6 = dL6_NN + blendW .* (dL6_Sig - dL6_NN);
end

function [dL4E, dL4I] = boundary_directional_finite_diff(celltype, L4EUse, L4IUse, L6Use, ContrastEff, OriCos2)
hE = 1e-4 .* max(1, abs(double(L4EUse)));
hI = 1e-4 .* max(1, abs(double(L4IUse)));
yEPlus = call_h96_pref6D_predictor(celltype, L4EUse + hE, L4IUse, L6Use, ContrastEff, OriCos2, false);
yEMinus = call_h96_pref6D_predictor(celltype, L4EUse - hE, L4IUse, L6Use, ContrastEff, OriCos2, false);
yIPlus = call_h96_pref6D_predictor(celltype, L4EUse, L4IUse + hI, L6Use, ContrastEff, OriCos2, false);
yIMinus = call_h96_pref6D_predictor(celltype, L4EUse, L4IUse - hI, L6Use, ContrastEff, OriCos2, false);
dL4E = (yEPlus - yEMinus) ./ (2 .* hE);
dL4I = (yIPlus - yIMinus) ./ (2 .* hI);
dL4E(~isfinite(dL4E)) = 0;
dL4I(~isfinite(dL4I)) = 0;
end

function [y, dL4E, dL4I, dL6] = call_h96_pref6D_predictor(celltype, L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, needGrad)
if needGrad
    switch celltype
        case 'S'
            [y, dL4E, dL4I, dL6] = predict_pref6D_S_with_grad(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2);
        case 'C'
            [y, dL4E, dL4I, dL6] = predict_pref6D_C_with_grad(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2);
        case 'I'
            [y, dL4E, dL4I, dL6] = predict_pref6D_I_with_grad(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2);
        otherwise
            error('Unknown celltype');
    end
else
    switch celltype
        case 'S'
            y = predict_pref6D_S(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2);
        case 'C'
            y = predict_pref6D_C(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2);
        case 'I'
            y = predict_pref6D_I(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2);
        otherwise
            error('Unknown celltype');
    end
    dL4E = zeros(size(y));
    dL4I = zeros(size(y));
    dL6 = zeros(size(y));
end
end
