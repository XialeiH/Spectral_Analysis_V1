function result = l6ns_inhibition_nearsilent_branch(setup,multistartResult,outputDir)
% Continue the near-silent inhibition root found at full inhibition freeze.

context = setup.Context;
fixed = context.FixedPoint(:);
phi = @(x,w)l6ns_phi(x,0,context,[1 1],1-w);
stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
weightScale = 0.1;
[~,rootIndex] = max(multistartResult.RootTable.distanceFromPersistent);
firstState = multistartResult.RootStates{rootIndex};
firstW = multistartResult.FreezeWeight;

initialStep = 0.01;
initialized = false;
while initialStep>=5e-4
    secondW = firstW-initialStep;
    [secondState,diagnostic] = l6ns_fixed_point_polish(phi,fixed,secondW,firstState);
    if diagnostic.converged
        initialized = true;
        break
    end
    initialStep = initialStep/2;
end
if ~initialized
    error('Could not initialize the near-silent branch below w=1.');
end

[tangentState,tangentW,lastDistance] = local_secant(firstState,firstW, ...
    secondState,secondW,stateScale,weightScale,[],[]);
states = {firstState;secondState};
weights = [firstW;secondW];
diagnostics = {struct('converged',true,'iterations',0, ...
    'residualNorm',norm(phi(firstState,firstW)-firstState), ...
    'relativeResidual',norm(phi(firstState,firstW)-firstState)/max(1,norm(fixed))); ...
    diagnostic};
step = min(0.5,max(0.05,0.8*lastDistance));
minimumStep = 0.003;
maximumStep = 0.5;
maximumPoints = 180;
attemptRows = {};
terminalReason = "maximum_points";
previousState = secondState;
previousW = secondW;

for pointIndex = 1:maximumPoints
    stateGuess = previousState+step*tangentState;
    wGuess = previousW+step*tangentW;
    [state,w,diagnostic] = local_arclength_newton(phi,previousState,previousW, ...
        tangentState,tangentW,step,stateGuess,wGuess,stateScale,weightScale,fixed);
    attemptRows(end+1,:) = {pointIndex,step,w,diagnostic.converged, ... %#ok<AGROW>
        diagnostic.residualNorm,diagnostic.relativeResidual,diagnostic.iterations, ...
        norm(state-fixed),min(state),max(state)};
    fprintf('Inhibition near-silent arc: point=%d ds=%.4g w=%+.9g dist=%.5g residual=%.3e converged=%d\n', ...
        pointIndex,step,w,norm(state-fixed),diagnostic.residualNorm,diagnostic.converged);
    if ~diagnostic.converged
        step = step/2;
        if step<minimumStep
            terminalReason = "minimum_step_failure";
            break
        end
        continue
    end

    states{end+1,1} = state; %#ok<AGROW>
    weights(end+1,1) = w; %#ok<AGROW>
    diagnostics{end+1,1} = diagnostic; %#ok<AGROW>
    [newTangentState,newTangentW,~] = local_secant(previousState,previousW, ...
        state,w,stateScale,weightScale,tangentState,tangentW);
    previousState = state;
    previousW = w;
    tangentState = newTangentState;
    tangentW = newTangentW;
    step = min(maximumStep,step*1.12);

    if w<-0.10
        terminalReason = "passed_w_minus_0p1";
        break
    elseif w>1.10
        terminalReason = "passed_w_1p1";
        break
    elseif pointIndex>12 && norm(state-fixed)<0.05*stateScale
        terminalReason = "returned_to_persistent_branch";
        break
    end
end

branchTable = local_branch_table(states,weights,diagnostics,fixed,context.MapSize);
samplePositions = unique(round(linspace(1,height(branchTable),min(13,height(branchTable)))));
stabilityRows = {};
for position = samplePositions
    modes = l6ns_numeric_leading_modes(@(x)phi(x,weights(position)), ...
        states{position},4,setup.Config);
    stabilityRows(end+1,:) = {position,weights(position), ... %#ok<AGROW>
        branchTable.relativeDistance(position),max(real(modes.Lambda))};
end
stabilityTable = cell2table(stabilityRows,'VariableNames', ...
    {'branchIndex','freezeWeight','relativeDistance','maxRealLambda'});
attemptTable = cell2table(attemptRows,'VariableNames',{'pointIndex','attemptedStep', ...
    'freezeWeight','converged','residualNorm','relativeResidual','newtonIterations', ...
    'distanceFromPersistent','minimumState','maximumState'});

if ~exist(outputDir,'dir'); mkdir(outputDir); end
writetable(branchTable,fullfile(outputDir,'inhibition_nearsilent_branch.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(stabilityTable,fullfile(outputDir,'inhibition_nearsilent_stability.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(attemptTable,fullfile(outputDir,'inhibition_nearsilent_attempts.tsv'), ...
    'FileType','text','Delimiter','\t');
result = struct('Table',branchTable,'States',{states},'Stability',stabilityTable, ...
    'Attempts',attemptTable,'TerminalReason',terminalReason);
save(fullfile(outputDir,'inhibition_nearsilent_branch_result.mat'),'result','-v7.3');
local_plot(result,setup.Inhibition.W,outputDir);
end

function tableOut = local_branch_table(states,weights,diagnostics,fixed,mapSize)
n = prod(mapSize);
rows = cell(numel(states),11);
for index = 1:numel(states)
    state = states{index};
    rows(index,:) = {index,weights(index),norm(state-fixed), ...
        norm(state-fixed)/norm(fixed),mean(state(1:n)),mean(state(n+1:2*n)), ...
        mean(state(2*n+1:3*n)),min(state),max(state), ...
        diagnostics{index}.residualNorm,diagnostics{index}.relativeResidual};
end
tableOut = cell2table(rows,'VariableNames',{'branchIndex','freezeWeight', ...
    'distanceFromPersistent','relativeDistance','meanS','meanC','meanI', ...
    'minimumState','maximumState','residualNorm','relativeResidual'});
end

function local_plot(result,criticalW,outputDir)
fig = figure('Visible','off','Color','w','Position',[80 50 1180 1000]);
layout = tiledlayout(3,1,'TileSpacing','loose','Padding','loose');
ax1 = nexttile(layout); plot(ax1,result.Table.freezeWeight, ...
    result.Table.relativeDistance,'o-','LineWidth',1.5,'MarkerSize',3);
xline(ax1,criticalW,'r:','LineWidth',1.2); grid(ax1,'on');
xlabel(ax1,'inhibition freeze weight w_I'); ylabel(ax1,'||f-f^*||_2 / ||f^*||_2');
title(ax1,'Near-silent fixed-point branch relative to the persistent operating point');
ax2 = nexttile(layout); plot(ax2,result.Table.freezeWeight, ...
    [result.Table.meanS result.Table.meanC result.Table.meanI],'-','LineWidth',1.5);
xline(ax2,criticalW,'r:','LineWidth',1.2); grid(ax2,'on');
xlabel(ax2,'inhibition freeze weight w_I'); ylabel(ax2,'mean firing rate');
legend(ax2,{'S','C','I'},'Location','best'); title(ax2,'Population means on near-silent branch');
ax3 = nexttile(layout); plot(ax3,result.Stability.freezeWeight, ...
    result.Stability.maxRealLambda,'o-','LineWidth',1.5);
yline(ax3,1,'r--','LineWidth',1.2); xline(ax3,criticalW,'r:','LineWidth',1.2);
grid(ax3,'on'); xlabel(ax3,'inhibition freeze weight w_I');
ylabel(ax3,'max Re \lambda(D\Phi)'); title(ax3,'Sampled stability of near-silent branch');
set([ax1 ax2 ax3],'FontSize',11);
l6ns_save_figure(fig,outputDir,'inhibition_nearsilent_global_branch'); close(fig);
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
