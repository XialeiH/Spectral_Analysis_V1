function result = l6ns_global_bifurcation_arclength(setup,refinedResult,outputDir,options)
% Follow a signed secondary branch through folds by pseudo-arclength continuation.

if nargin<4; options = struct(); end

name = refinedResult.Name;
branchSign = refinedResult.Sign;
critical = setup.(name);
context = setup.Context;
fixed = context.FixedPoint(:);
if strcmp(name,'L6')
    phi = @(x,w)l6ns_phi(x,w,context,[1 1],1);
else
    phi = @(x,w)l6ns_phi(x,0,context,[1 1],1-w);
end

valid = find(refinedResult.Table.converged);
if numel(valid)<2
    error('At least two converged refined points are required.');
end
tableOut = refinedResult.Table(valid,:);
states = refinedResult.States(valid);
stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
weightScale = 0.1;

olderState = states{end-1};
olderW = tableOut.freezeWeight(end-1);
previousState = states{end};
previousW = tableOut.freezeWeight(end);
[tangentState,tangentW,lastDistance] = local_secant(olderState,olderW, ...
    previousState,previousW,stateScale,weightScale,[],[]);
minimumStep = local_option(options,'MinimumStep',0.004);
maximumStep = local_option(options,'MaximumStep',0.8);
maximumPoints = local_option(options,'MaximumPoints',180);
returnDistanceRatio = local_option(options,'ReturnDistanceRatio',0.05);
step = min(maximumStep,max(min(0.15,maximumStep),0.6*lastDistance));
attemptRows = {};
terminalReason = "maximum_points";

for pointIndex = 1:maximumPoints
    stateGuess = previousState+step*tangentState;
    wGuess = previousW+step*tangentW;
    [state,w,diagnostic] = local_arclength_newton(phi,previousState,previousW, ...
        tangentState,tangentW,step,stateGuess,wGuess,stateScale,weightScale,fixed);
    amplitude = real(critical.Left'*(state-fixed));
    attemptRows(end+1,:) = {pointIndex,step,w,amplitude/stateScale, ... %#ok<AGROW>
        diagnostic.converged,diagnostic.residualNorm,diagnostic.relativeResidual, ...
        diagnostic.iterations,min(state),max(state)};
    fprintf('%s sign %+d arc: point=%d ds=%.4g w=%+.9g a/scale=%+.5g residual=%.3e converged=%d\n', ...
        name,branchSign,pointIndex,step,w,amplitude/stateScale, ...
        diagnostic.residualNorm,diagnostic.converged);
    if ~diagnostic.converged
        step = step/2;
        if step<minimumStep
            terminalReason = "minimum_step_failure";
            break
        end
        continue
    end

    newRow = table(amplitude,amplitude/stateScale,w,norm(state-fixed), ...
        diagnostic.residualNorm,diagnostic.relativeResidual,true, ...
        diagnostic.iterations,min(state),max(state), ...
        'VariableNames',tableOut.Properties.VariableNames);
    tableOut = [tableOut;newRow]; %#ok<AGROW>
    states{end+1,1} = state; %#ok<AGROW>
    [newTangentState,newTangentW,~] = local_secant(previousState,previousW, ...
        state,w,stateScale,weightScale,tangentState,tangentW);
    olderState = previousState; %#ok<NASGU>
    olderW = previousW; %#ok<NASGU>
    previousState = state;
    previousW = w;
    tangentState = newTangentState;
    tangentW = newTangentW;
    step = min(maximumStep,step*1.12);

    if min(state)<0
        terminalReason = "left_nonnegative_state_region";
        break
    elseif w>1.10
        terminalReason = "passed_w_1p1";
        break
    elseif w<-0.50
        terminalReason = "passed_w_minus_0p5";
        break
    elseif pointIndex>12 && norm(state-fixed)<returnDistanceRatio*stateScale
        terminalReason = "returned_to_persistent_branch";
        break
    end
end

attemptTable = cell2table(attemptRows,'VariableNames',{'pointIndex', ...
    'attemptedStep','freezeWeight','normalizedAmplitude','converged', ...
    'residualNorm','relativeResidual','newtonIterations','minimumState','maximumState'});
result = refinedResult;
result.Table = tableOut;
result.States = states;
result.ArclengthAttempts = attemptTable;
result.TerminalReason = terminalReason;

if ~exist(outputDir,'dir'); mkdir(outputDir); end
stem = sprintf('%s_sign_%+d',lower(name),branchSign);
stem = strrep(stem,'+','p');
stem = strrep(stem,'-','m');
writetable(tableOut,fullfile(outputDir,[stem '_continuation.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(attemptTable,fullfile(outputDir,[stem '_arclength_attempts.tsv']), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputDir,[stem '_result.mat']),'result','-v7.3');
end

function value = local_option(options,name,defaultValue)
if isfield(options,name)
    value = options.(name);
else
    value = defaultValue;
end
end

function [tangentState,tangentW,distance] = local_secant( ...
        firstState,firstW,secondState,secondW,stateScale,weightScale,oldState,oldW)
deltaState = secondState-firstState;
deltaW = secondW-firstW;
distance = sqrt(norm(deltaState/stateScale)^2+(deltaW/weightScale)^2);
tangentState = deltaState/distance;
tangentW = deltaW/distance;
if ~isempty(oldState)
    alignment = real(oldState'*tangentState)/(stateScale^2)+ ...
        oldW*tangentW/(weightScale^2);
    if alignment<0
        tangentState = -tangentState;
        tangentW = -tangentW;
    end
end
end

function [state,w,diagnostic] = local_arclength_newton(phi,baseState,baseW, ...
        tangentState,tangentW,targetDistance,state,w,stateScale,weightScale,fixed)
tolerance = 2e-9*max(1,norm(fixed));
residualNorm = inf;
converged = false;
iteration = 0;
for iteration = 1:30
    phiBase = phi(state,w);
    arcResidual = real(tangentState'*(state-baseState))/(stateScale^2)+ ...
        tangentW*(w-baseW)/(weightScale^2)-targetDistance;
    residual = [phiBase-state;arcResidual];
    residualNorm = norm(residual);
    if residualNorm<tolerance
        converged = true;
        break
    end
    stateStep = 4e-7*stateScale;
    parameterStep = 4e-7;
    dGdw = (phi(state,w+parameterStep)-phi(state,w-parameterStep))/(2*parameterStep);
    operator = @(delta)local_arclength_action(delta,state,w,phi,phiBase, ...
        stateStep,dGdw,tangentState,tangentW,stateScale,weightScale);
    [delta,~] = gmres(operator,-residual,[],2e-6,90);
    if any(~isfinite(delta)); break; end
    accepted = false;
    lineScale = 1;
    for lineIteration = 1:12
        candidateState = state+lineScale*delta(1:end-1);
        candidateW = w+lineScale*delta(end);
        candidateArc = real(tangentState'*(candidateState-baseState))/(stateScale^2)+ ...
            tangentW*(candidateW-baseW)/(weightScale^2)-targetDistance;
        candidateResidual = [phi(candidateState,candidateW)-candidateState;candidateArc];
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

function output = local_arclength_action(delta,state,w,phi,phiBase,stateStep, ...
        dGdw,tangentState,tangentW,stateScale,weightScale)
direction = delta(1:end-1);
parameterDirection = delta(end);
directionNorm = norm(direction);
if directionNorm==0
    dGdf = zeros(size(direction));
else
    h = stateStep/directionNorm;
    dGdf = (phi(state+h*direction,w)-phiBase)/h-direction;
end
arcDerivative = real(tangentState'*direction)/(stateScale^2)+ ...
    tangentW*parameterDirection/(weightScale^2);
output = [dGdf+dGdw*parameterDirection;arcDerivative];
end
