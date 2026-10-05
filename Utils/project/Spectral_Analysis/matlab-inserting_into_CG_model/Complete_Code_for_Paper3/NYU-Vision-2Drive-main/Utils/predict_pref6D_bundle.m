function [rate_S, rate_C, rate_I] = predict_pref6D_bundle(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, OriSin2) %#codegen
rate_S = predict_pref6D_S(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, OriSin2);
rate_C = predict_pref6D_C(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, OriSin2);
rate_I = predict_pref6D_I(L4EUse, L4IUse, L6Use, ContrastEff, OriCos2, OriSin2);
end
