function [J,operators,rowSumError] = frozen_kernel_jacobian(frozen,baseline,l4Alpha,l6Alpha)
% Assemble with saved local derivatives; no response function is evaluated here.
n = numel(frozen.FixedPoint)/3;
names = {'C_SS','C_SC','C_SI','C_CS','C_CC','C_CI','C_IS','C_IC','C_II','L6'};
operators = baseline;
rowSumError = zeros(numel(names),1);
for index = 1:numel(names)
    name = names{index};
    matrix = baseline.(name);
    rowSums = full(sum(matrix,2));
    assert(max(abs(rowSums-mean(rowSums))) <= 1e-10*max(1,max(abs(rowSums))), ...
        'Expected constant row sums for a spatial convolution operator.');
    alpha = l4Alpha;
    if strcmp(name,'L6'); alpha = l6Alpha; end
    operators.(name) = (1-alpha)*matrix+alpha*mean(rowSums)*speye(n);
    rowSumError(index) = max(abs(full(sum(operators.(name),2))-rowSums)) ...
        /max(max(abs(rowSums)),eps);
end
assert(max(rowSumError) < 1e-12, 'Kernel replacement must preserve row sums.');

d = frozen.Preprocess;
p = frozen.Normalization;
g = frozen.Response;
cSS = operators.C_SS; cSC = operators.C_SC; cSI = operators.C_SI;
cCS = operators.C_CS; cCC = operators.C_CC; cCI = operators.C_CI;
cIS = operators.C_IS; cIC = operators.C_IC; cII = operators.C_II;
dL4ESS = (cSS*d.SS+cSC*d.CS)/p.L4SEp;
dL4ESC = (cSS*d.SC+cSC*d.CC)/p.L4SEp;
dL4ECS = (cCS*d.SS+cCC*d.CS)/p.L4CEp;
dL4ECC = (cCS*d.SC+cCC*d.CC)/p.L4CEp;
dL4EIS = (cIS*d.SS+cIC*d.CS)/p.L4IEp;
dL4EIC = (cIS*d.SC+cIC*d.CC)/p.L4IEp;
dL4ISI = (cSI*d.II)/p.L4SIp;
dL4ICI = (cCI*d.II)/p.L4CIp;
dL4III = (cII*d.II)/p.L4IIp;
wS = 1-frozen.CWeight; wC = frozen.CWeight;
dEuseS = wS*d.SS+wC*d.CS;
dEuseC = wS*d.SC+wC*d.CC;
dL6S = frozen.L6Derivative*operators.L6*dEuseS;
dL6C = frozen.L6Derivative*operators.L6*dEuseC;
J = [g.S4E*dL4ESS+g.S6*dL6S, g.S4E*dL4ESC+g.S6*dL6C, g.S4I*dL4ISI; ...
     g.C4E*dL4ECS+g.C6*dL6S, g.C4E*dL4ECC+g.C6*dL6C, g.C4I*dL4ICI; ...
     g.I4E*dL4EIS+g.I6*dL6S, g.I4E*dL4EIC+g.I6*dL6C, g.I4I*dL4III];
end
