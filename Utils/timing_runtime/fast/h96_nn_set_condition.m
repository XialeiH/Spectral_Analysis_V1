function fast = h96_nn_set_condition(base, contrast, orientation)
% Add only the compact static predictor inputs for one stimulus condition.
n = base.n;
alpha = orientation(:);
if isscalar(alpha)
    alpha = repmat(alpha, n, 1);
end
contrastVector = contrast(:);
if isscalar(contrastVector)
    contrastVector = repmat(contrastVector, n, 1);
end
assert(numel(alpha)==n && numel(contrastVector)==n, ...
    'Contrast and orientation must be scalar or N-by-1.');

rows = base.mixRows;
categories = base.mixCategories;
gamma = mod(alpha(rows)-base.thetaPref(categories).'+90, 180)-90;
orientationFeature = cosd(2.*abs(gamma));
contrastFeature = contrastVector(rows);
contrastFeature(categories==5) = 0;

fast = base;
fast.staticCompactPages = zeros(base.activePairCount, 5, 3, base.nnPrecision);
for population = 1:3
    fast.staticCompactPages(:,1,population) = ...
        (cast(orientationFeature, base.nnPrecision)-base.mu(1,1,population))./base.sd(1,1,population);
    fast.staticCompactPages(:,2,population) = ...
        (cast(contrastFeature, base.nnPrecision)-base.mu(1,2,population))./base.sd(1,2,population);
end
end
