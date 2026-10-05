function [state,diagnostic] = l6ns_fixed_point_newton(phi,fixed,w,state)
% Matrix-free damped Newton solve of Phi(x,w)-x=0.

stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
tolerance = 2e-9*max(1,norm(fixed));
residualNorm = inf;
converged = false;
iteration = 0;
for iteration = 1:24
    phiBase = phi(state,w);
    residual = phiBase-state;
    residualNorm = norm(residual);
    if residualNorm<tolerance
        converged = true;
        break
    end
    stateStep = 4e-7*stateScale;
    operator = @(direction)local_action(direction,state,w,phi,phiBase,stateStep);
    [delta,~] = gmres(operator,-residual,[],2e-6,80);
    if any(~isfinite(delta)); break; end
    accepted = false;
    lineScale = 1;
    for lineIteration = 1:10
        candidate = state+lineScale*delta;
        if norm(phi(candidate,w)-candidate)<residualNorm
            state = candidate;
            accepted = true;
            break
        end
        lineScale = lineScale/2;
    end
    if ~accepted; break; end
end
diagnostic = struct('converged',converged,'iterations',iteration, ...
    'residualNorm',residualNorm, ...
    'relativeResidual',residualNorm/max(1,norm(fixed)));
end

function output = local_action(direction,state,w,phi,phiBase,stateStep)
directionNorm = norm(direction);
if directionNorm==0
    output = zeros(size(direction));
    return
end
h = stateStep/directionNorm;
output = (phi(state+h*direction,w)-phiBase)/h-direction;
end
