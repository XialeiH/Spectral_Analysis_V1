function result = l6ns_global_bifurcation_refine(setup,coarseResult,outputDir)
% Continue a coarse signed branch with adaptive amplitude steps.

name = coarseResult.Name;
branchSign = coarseResult.Sign;
critical = setup.(name);
context = setup.Context;
fixed = context.FixedPoint(:);
if strcmp(name,'L6')
    phi = @(x,w)l6ns_phi(x,w,context,[1 1],1);
else
    phi = @(x,w)l6ns_phi(x,0,context,[1 1],1-w);
end

valid = find(coarseResult.Table.converged);
if numel(valid)<2
    error('At least two converged coarse points are required.');
end
tableOut = coarseResult.Table(valid,:);
states = coarseResult.States(valid);
stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));

olderTarget = tableOut.amplitude(end-1);
olderState = states{end-1};
olderW = tableOut.freezeWeight(end-1);
previousTarget = tableOut.amplitude(end);
previousState = states{end};
previousW = tableOut.freezeWeight(end);

step = 0.35;
minimumStep = 0.01;
maximumStep = 0.75;
maximumAmplitude = 60;
attemptRows = {};
terminalReason = "maximum_amplitude";
attempt = 0;
while abs(previousTarget/stateScale)<maximumAmplitude
    attempt = attempt+1;
    nextNormalized = min(maximumAmplitude,abs(previousTarget/stateScale)+step);
    target = branchSign*stateScale*nextNormalized;
    ratio = (target-previousTarget)/(previousTarget-olderTarget);
    stateGuess = previousState+ratio*(previousState-olderState);
    wGuess = previousW+ratio*(previousW-olderW);
    [state,w,diagnostic] = local_augmented_newton(phi,fixed,critical.Left, ...
        target,stateGuess,wGuess,stateScale);
    attemptRows(end+1,:) = {attempt,target/stateScale,step,w, ... %#ok<AGROW>
        diagnostic.converged,diagnostic.residualNorm,diagnostic.relativeResidual, ...
        diagnostic.iterations,min(state),max(state)};
    fprintf('%s sign %+d refine: a/scale=%+.5g step=%.4g w=%+.9g residual=%.3e converged=%d\n', ...
        name,branchSign,target/stateScale,step,w,diagnostic.residualNorm,diagnostic.converged);
    if ~diagnostic.converged
        step = step/2;
        if step<minimumStep
            terminalReason = "minimum_step_failure";
            break
        end
        continue
    end

    newRow = table(target,target/stateScale,w,norm(state-fixed), ...
        diagnostic.residualNorm,diagnostic.relativeResidual,true, ...
        diagnostic.iterations,min(state),max(state), ...
        'VariableNames',tableOut.Properties.VariableNames);
    tableOut = [tableOut;newRow]; %#ok<AGROW>
    states{end+1,1} = state; %#ok<AGROW>
    olderTarget = previousTarget;
    olderState = previousState;
    olderW = previousW;
    previousTarget = target;
    previousState = state;
    previousW = w;
    step = min(maximumStep,step*1.15);

    if min(state)<0
        terminalReason = "left_nonnegative_state_region";
        break
    elseif strcmp(name,'L6') && w>1.10
        terminalReason = "passed_w_1p1";
        break
    elseif strcmp(name,'Inhibition') && w<-0.10
        terminalReason = "passed_w_minus_0p1";
        break
    elseif nextNormalized>=maximumAmplitude
        terminalReason = "maximum_amplitude";
        break
    end
end

attemptTable = cell2table(attemptRows,'VariableNames',{'attempt', ...
    'normalizedAmplitude','attemptedStep','freezeWeight','converged', ...
    'residualNorm','relativeResidual','newtonIterations','minimumState','maximumState'});
result = coarseResult;
result.Table = tableOut;
result.States = states;
result.RefinementAttempts = attemptTable;
result.TerminalReason = terminalReason;

if ~exist(outputDir,'dir'); mkdir(outputDir); end
stem = sprintf('%s_sign_%+d',lower(name),branchSign);
stem = strrep(stem,'+','p');
stem = strrep(stem,'-','m');
writetable(tableOut,fullfile(outputDir,[stem '_continuation.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(attemptTable,fullfile(outputDir,[stem '_refinement_attempts.tsv']), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputDir,[stem '_result.mat']),'result','-v7.3');
end

function [state,w,diagnostic] = local_augmented_newton( ...
        phi,fixed,left,target,state,w,stateScale)
tolerance = 2e-9*max(1,norm(fixed));
residualNorm = inf;
converged = false;
iteration = 0;
for iteration = 1:30
    phiBase = phi(state,w);
    residual = [phiBase-state;left'*(state-fixed)-target];
    residualNorm = norm(residual);
    if residualNorm<tolerance
        converged = true;
        break
    end
    stateStep = 4e-7*stateScale;
    parameterStep = 4e-7;
    dGdw = (phi(state,w+parameterStep)-phi(state,w-parameterStep))/(2*parameterStep);
    operator = @(delta)local_augmented_action(delta,state,w,phi,phiBase, ...
        stateStep,dGdw,left);
    [delta,~] = gmres(operator,-residual,[],2e-6,90);
    if any(~isfinite(delta)); break; end
    accepted = false;
    lineScale = 1;
    for lineIteration = 1:12
        candidateState = state+lineScale*delta(1:end-1);
        candidateW = w+lineScale*delta(end);
        candidateResidual = [phi(candidateState,candidateW)-candidateState; ...
            left'*(candidateState-fixed)-target];
        if norm(candidateResidual)<residualNorm
            state = candidateState;
            w = candidateW;
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

function output = local_augmented_action(delta,state,w,phi,phiBase,stateStep,dGdw,left)
direction = delta(1:end-1);
parameterDirection = delta(end);
directionNorm = norm(direction);
if directionNorm==0
    dGdf = zeros(size(direction));
else
    h = stateStep/directionNorm;
    dGdf = (phi(state+h*direction,w)-phiBase)/h-direction;
end
output = [dGdf+dGdw*parameterDirection;left'*direction];
end
