function [Sigma, rhat] = GA_empirical_covariance(Y, shrinkage)
% Sample covariance with optional diagonal shrinkage from proposal section 5.4.
% Rows are response variables; columns are repeated trials.
if nargin < 2 || isempty(shrinkage)
    shrinkage = 0;
end

rhat = mean(Y, 2);
Y0 = Y - rhat;
nTrials = size(Y, 2);
SigmaRaw = (Y0 * Y0') ./ max(nTrials - 1, 1);
diagTarget = trace(SigmaRaw) / size(SigmaRaw, 1) * eye(size(SigmaRaw, 1));
Sigma = (1 - shrinkage) * SigmaRaw + shrinkage * diagTarget;
Sigma = 0.5 * (Sigma + Sigma');
end
