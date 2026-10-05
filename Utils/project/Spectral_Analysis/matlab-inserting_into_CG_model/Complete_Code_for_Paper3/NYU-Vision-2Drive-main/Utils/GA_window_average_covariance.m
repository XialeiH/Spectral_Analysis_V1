function SigmaY = GA_window_average_covariance(A, SigmaState, nSteps)
% Discrete covariance of a window-averaged observable, proposal section 5.7.
% Assumes stationary linearized dynamics with lag covariance Gamma(k)=A^k*Sigma.
if nargin < 3 || isempty(nSteps)
    nSteps = 1;
end

SigmaY = nSteps * SigmaState;
Ak = speye(size(A, 1));
for k = 1:nSteps-1
    Ak = A * Ak;
    gammaK = Ak * SigmaState;
    SigmaY = SigmaY + (nSteps - k) * (gammaK + gammaK');
end
SigmaY = SigmaY ./ (nSteps^2);
SigmaY = 0.5 * (SigmaY + SigmaY');
end
