function result = l6ns_global_bifurcation_branch(setup,taskId,outputDir)
% Continue one signed secondary branch well beyond the local normal form.

names = {'L6','L6','Inhibition','Inhibition'};
signs = [-1 1 -1 1];
if taskId<1 || taskId>numel(names)
    error('taskId must be in 1:4.');
end
name = names{taskId};
branchSign = signs(taskId);
critical = setup.(name);
context = setup.Context;
fixed = context.FixedPoint(:);
cfg = setup.Config;
if strcmp(name,'L6')
    phi = @(x,w)l6ns_phi(x,w,context,[1 1],1);
else
    phi = @(x,w)l6ns_phi(x,0,context,[1 1],1-w);
end

stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
normalized = [0.004 0.008 0.015 0.025];
value = 0.04;
while value<=80
    normalized(end+1) = value; %#ok<AGROW>
    value = value*1.28;
end
targets = branchSign*stateScale*normalized;

rows = {};
states = {};
previousTarget = 0;
previousState = fixed;
previousW = critical.W;
olderTarget = NaN;
olderState = [];
olderW = NaN;
for targetIndex = 1:numel(targets)
    target = targets(targetIndex);
    if isempty(olderState)
        state = fixed+target*critical.Right;
        w = critical.W-real((critical.Quadratic*target+ ...
            critical.Cubic*target^2)/critical.Beta);
    else
        ratio = (target-previousTarget)/(previousTarget-olderTarget);
        state = previousState+ratio*(previousState-olderState);
        w = previousW+ratio*(previousW-olderW);
    end
    [state,w,diagnostic] = local_augmented_newton(phi,fixed,critical.Left, ...
        target,state,w,stateScale);
    rows(end+1,:) = {target,target/stateScale,w,norm(state-fixed), ...
        diagnostic.residualNorm,diagnostic.relativeResidual, ...
        diagnostic.converged,diagnostic.iterations,min(state),max(state)}; %#ok<AGROW>
    states{end+1,1} = state; %#ok<AGROW>
    fprintf('%s sign %+d: a/scale=%+.5g w=%+.9g residual=%.3e converged=%d\n', ...
        name,branchSign,target/stateScale,w,diagnostic.residualNorm,diagnostic.converged);
    if ~diagnostic.converged
        break
    end
    olderTarget = previousTarget;
    olderState = previousState;
    olderW = previousW;
    previousTarget = target;
    previousState = state;
    previousW = w;
    if strcmp(name,'L6') && w>1.10
        break
    elseif strcmp(name,'Inhibition') && w<-0.10
        break
    end
end

tableOut = cell2table(rows,'VariableNames',{'amplitude','normalizedAmplitude', ...
    'freezeWeight','branchDistance','residualNorm','relativeResidual', ...
    'converged','newtonIterations','minimumState','maximumState'});
targetRoots = local_target_roots(phi,fixed,tableOut,states,cfg,[0 1]);

if ~exist(outputDir,'dir'); mkdir(outputDir); end
stem = sprintf('%s_sign_%+d',lower(name),branchSign);
stem = strrep(stem,'+','p');
stem = strrep(stem,'-','m');
writetable(tableOut,fullfile(outputDir,[stem '_continuation.tsv']), ...
    'FileType','text','Delimiter','\t');
writetable(targetRoots.Table,fullfile(outputDir,[stem '_target_roots.tsv']), ...
    'FileType','text','Delimiter','\t');
result = struct('Name',name,'Sign',branchSign,'Critical',critical, ...
    'Table',tableOut,'States',{states},'TargetRoots',targetRoots);
save(fullfile(outputDir,[stem '_result.mat']),'result','-v7.3');
end

function [state,w,diagnostic] = local_augmented_newton( ...
        phi,fixed,left,target,state,w,stateScale)
tolerance = 2e-9*max(1,norm(fixed));
residualNorm = inf;
converged = false;
iteration = 0;
for iteration = 1:22
    phiBase = phi(state,w);
    residual = [phiBase-state;left'*(state-fixed)-target];
    residualNorm = norm(residual);
    if residualNorm<tolerance
        converged = true;
        break
    end
    stateStep = 4e-7*stateScale;
    parameterStep = 4e-7;
    dGdw = (phi(state,w+parameterStep)-phi(state,w-parameterStep)) / ...
        (2*parameterStep);
    operator = @(delta)local_augmented_action(delta,state,w,phi,phiBase, ...
        stateStep,dGdw,left);
    [delta,flag] = gmres(operator,-residual,[],2e-6,70); %#ok<ASGLU>
    if any(~isfinite(delta)); break; end
    accepted = false;
    lineScale = 1;
    for lineIteration = 1:10
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

function roots = local_target_roots(phi,fixed,tableOut,states,cfg,targetWeights)
rows = {};
rootStates = {};
valid = find(tableOut.converged);
for targetIndex = 1:numel(targetWeights)
    targetW = targetWeights(targetIndex);
    bracket = [];
    for index = 1:numel(valid)-1
        first = valid(index);
        second = valid(index+1);
        if (tableOut.freezeWeight(first)-targetW)* ...
                (tableOut.freezeWeight(second)-targetW)<=0
            bracket = [first second]; %#ok<AGROW>
            break
        end
    end
    if isempty(bracket)
        rows(end+1,:) = {targetW,false,NaN,NaN,NaN,NaN,NaN}; %#ok<AGROW>
        rootStates{end+1,1} = []; %#ok<AGROW>
        continue
    end
    w1 = tableOut.freezeWeight(bracket(1));
    w2 = tableOut.freezeWeight(bracket(2));
    fraction = (targetW-w1)/(w2-w1);
    state = states{bracket(1)}+fraction*(states{bracket(2)}-states{bracket(1)});
    [state,diagnostic] = l6ns_fixed_point_newton(phi,fixed,targetW,state);
    maxReal = NaN;
    if diagnostic.converged
        modes = l6ns_numeric_leading_modes(@(x)phi(x,targetW),state,4,cfg);
        maxReal = max(real(modes.Lambda));
    end
    rows(end+1,:) = {targetW,diagnostic.converged,norm(state-fixed), ...
        diagnostic.residualNorm,diagnostic.relativeResidual,maxReal,min(state)}; %#ok<AGROW>
    rootStates{end+1,1} = state; %#ok<AGROW>
end
roots = struct();
roots.Table = cell2table(rows,'VariableNames',{'freezeWeight','converged', ...
    'branchDistance','residualNorm','relativeResidual','maxRealLambda','minimumState'});
roots.States = rootStates;
end
