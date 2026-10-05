function k = h96_prepare_ekp_fast(EKp)
% Compile the active quadratic E saturation conversion.
assert(strcmp(EKp{end}, 'quadratic'), ...
    'Fast EKp currently requires the active quadratic mode.');
xData = EKp{4};
yData = EKp{3};
A = [xData.^2, xData, ones(3, 1)];
k.Coeff = A\yData;
k.Threshold = yData(1);
end
