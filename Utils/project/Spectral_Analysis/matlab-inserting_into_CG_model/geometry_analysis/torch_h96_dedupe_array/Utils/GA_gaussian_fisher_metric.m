function F = GA_gaussian_fisher_metric(Dr, Sigma, dSigma, reg)
% Full Gaussian Fisher correction from proposal section 5.10.
if nargin < 4 || isempty(reg)
    reg = 1e-6;
end

Prec = GA_regularized_precision(Sigma, reg);
nCoord = size(Dr, 2);
F = Dr' * Prec * Dr;

if nargin >= 3 && ~isempty(dSigma)
    for i = 1:nCoord
        for j = i:nCoord
            corr = 0.5 * trace(Prec * dSigma{i} * Prec * dSigma{j});
            F(i,j) = F(i,j) + corr;
            F(j,i) = F(i,j);
        end
    end
end
end
