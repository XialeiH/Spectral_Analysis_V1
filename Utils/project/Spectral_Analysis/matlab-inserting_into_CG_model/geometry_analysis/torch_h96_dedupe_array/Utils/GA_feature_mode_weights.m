function alpha = GA_feature_mode_weights(V, W, lambda, dPhi)
% Mode weights alpha_{k,i} for feature derivatives dPhi(:,i).
nModes = numel(lambda);
nFeatures = size(dPhi, 2);
alpha = zeros(nModes, nFeatures);
gap = abs(1 - lambda(:)).^2;
gap(gap < eps) = eps;

for i = 1:nFeatures
    drive = abs(W' * dPhi(:,i)).^2 ./ gap;
    total = sum(drive);
    if total > 0
        alpha(:,i) = drive ./ total;
    end
end
end
