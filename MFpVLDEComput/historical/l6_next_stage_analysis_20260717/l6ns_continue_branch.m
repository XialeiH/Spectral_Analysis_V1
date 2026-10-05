function continuation = l6ns_continue_branch(context, critical, cfg)
% Amplitude-constrained Newton-GMRES continuation of the second branch.

fixed = context.FixedPoint;
r = real(critical.Right);
l = real(critical.Left);
r = r / norm(r);
l = l / (l' * r);
scale = max(1, norm(fixed)/sqrt(numel(fixed)));
amplitudes = scale * [-0.025 -0.015 -0.008 -0.004 0.004 0.008 0.015 0.025];

rows = cell(numel(amplitudes), 1);
states = cell(numel(amplitudes), 1);
for ai = 1:numel(amplitudes)
    targetAmplitude = amplitudes(ai);
    w = critical.W - real(critical.Quadratic/critical.Beta) * targetAmplitude;
    state = fixed + targetAmplitude * r;
    converged = false;
    residualNorm = inf;
    convergenceTolerance = 5e-10 * max(1,norm(fixed));
    for newtonIteration = 1:18
        phi = @(x,weight) l6ns_phi(x,weight,context);
        g = phi(state,w) - state;
        residual = [g; l'*(state-fixed)-targetAmplitude];
        residualNorm = norm(residual);
        if residualNorm < convergenceTolerance
            converged = true;
            break
        end
        stateStep = 5e-7 * max(1,norm(state)/sqrt(numel(state)));
        weightStep = 5e-7;
        phiBase = phi(state,w);
        dGdw = (phi(state,w+weightStep)-phi(state,w-weightStep))/(2*weightStep);
        operator = @(delta) local_augmented_action(delta,state,w,phi,phiBase, ...
            stateStep,dGdw,l);
        [delta,flag] = gmres(operator,-residual,[],1e-6,60); %#ok<ASGLU>
        if any(~isfinite(delta)); break; end
        lineScale = 1;
        accepted = false;
        for lineIteration = 1:8
            candidateState = state + lineScale*delta(1:end-1);
            candidateW = w + lineScale*delta(end);
            candidateResidual = [phi(candidateState,candidateW)-candidateState; ...
                l'*(candidateState-fixed)-targetAmplitude];
            if norm(candidateResidual) < residualNorm
                state = candidateState;
                w = candidateW;
                accepted = true;
                break
            end
            lineScale = lineScale/2;
        end
        if ~accepted; break; end
    end
    maxReal = NaN;
    if converged
        phiAtW = @(x) l6ns_phi(x,w,context);
        nonlinearModes = l6ns_numeric_leading_modes(phiAtW,state,4,cfg);
        maxReal = max(real(nonlinearModes.Lambda));
    end
    rows{ai} = {targetAmplitude,w,norm(state-fixed),residualNorm, ...
        residualNorm/max(1,norm(fixed)),converged,newtonIteration,maxReal};
    states{ai} = state;
    fprintf('Continuation a=%+.5g w=%+.6g residual %.3e converged=%d.\n', ...
        targetAmplitude,w,residualNorm,converged);
end

continuation = struct();
continuation.Table = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'amplitude','w','branchDistance','residualNorm','relativeResidual', ...
    'converged','newtonIterations','maxRealLambda'});
continuation.States = states;
continuation.FixedPoint = fixed;
end

function output = local_augmented_action(delta,state,w,phi,phiBase,stateStep,dGdw,l)
direction = delta(1:end-1);
weightDirection = delta(end);
directionNorm = norm(direction);
if directionNorm == 0
    dGdf = zeros(size(direction));
else
    h = stateStep / directionNorm;
    dPhi = (phi(state+h*direction,w)-phiBase)/h;
    dGdf = dPhi-direction;
end
output = [dGdf+dGdw*weightDirection; l'*direction];
end
