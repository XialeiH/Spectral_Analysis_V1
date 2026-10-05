function Sigma = GA_lowrank_covariance(diagVar, U, lambda, noiseFloor)
% Low-rank-plus-diagonal covariance from proposal section 5.5.
if nargin < 4 || isempty(noiseFloor)
    noiseFloor = 0;
end

diagVar = diagVar(:) + noiseFloor;
D = spdiags(diagVar, 0, numel(diagVar), numel(diagVar));

if isempty(U)
    Sigma = D;
    return
end

if isvector(lambda)
    Lambda = diag(lambda(:));
else
    Lambda = lambda;
end

Sigma = D + U * Lambda * U';
Sigma = 0.5 * (Sigma + Sigma');
end
