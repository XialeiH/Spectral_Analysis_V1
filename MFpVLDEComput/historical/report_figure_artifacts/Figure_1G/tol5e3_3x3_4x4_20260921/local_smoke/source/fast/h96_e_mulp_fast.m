function y = h96_e_mulp_fast(x, k)
% Active quadratic E saturation conversion with precomputed coefficients.
y = x;
index = x>k.Threshold;
a = k.Coeff(1);
b = k.Coeff(2);
c = k.Coeff(3);
y(index) = (-b+sqrt(b^2-4*a*c+4*a*x(index)))/(2*a);
end
