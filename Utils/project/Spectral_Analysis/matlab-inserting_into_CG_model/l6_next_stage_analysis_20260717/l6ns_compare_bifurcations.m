function result = l6ns_compare_bifurcations(cfg,data,context,l6Result,primary)
% Compare the local nonlinear bifurcations produced by L6 and inhibition.

outputDir = fullfile(cfg.OutputRoot,'l6_inhibition_bifurcation');
if ~exist(outputDir,'dir'); mkdir(outputDir); end

pathway = primary.Pathway;
fixed = context.FixedPoint(:);
rawContext = context;
rawContext.FixedPointCorrection = [];
context.FixedPointCorrection = fixed-l6ns_phi(fixed,0,rawContext,[1 1],1);

l6 = local_existing_l6_result(cfg,pathway,l6Result);
inhibition = local_inhibition_result(cfg,pathway,context,fixed);

w6Grid = linspace(-0.35,0.05,81)';
wIGrid = linspace(-0.05,0.30,71)';
l6.StabilityCurve = local_stability_curve(cfg,pathway,w6Grid,'L6');
inhibition.StabilityCurve = local_stability_curve(cfg,pathway,wIGrid,'I');

writetable(l6.Continuation,fullfile(outputDir,'l6_secondary_branch.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(inhibition.Continuation,fullfile(outputDir,'inhibition_secondary_branch.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(l6.StabilityCurve,fullfile(outputDir,'l6_persistent_branch_stability.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(inhibition.StabilityCurve,fullfile(outputDir,'inhibition_persistent_branch_stability.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(inhibition.NormalForm,fullfile(outputDir,'inhibition_normal_form_convergence.tsv'), ...
    'FileType','text','Delimiter','\t');

summary = table(["L6";"inhibition"], ...
    [l6.CriticalFreezeWeight;inhibition.CriticalFreezeWeight], ...
    [l6.CriticalGain;inhibition.CriticalGain], ...
    [l6.Quadratic;inhibition.Quadratic], ...
    [l6.Cubic;inhibition.Cubic], ...
    [l6.Beta;inhibition.Beta], ...
    [l6.ZeroResidual;inhibition.ZeroResidual], ...
    [l6.StabilityExchange;inhibition.StabilityExchange], ...
    [string(l6.Classification);string(inhibition.Classification)], ...
    'VariableNames',{'pathway','criticalFreezeWeight','criticalDynamicGain', ...
    'quadraticCoefficient','cubicCoefficient','parameterSensitivity', ...
    'zeroStateResidual','stabilityExchange','classification'});
writetable(summary,fullfile(outputDir,'bifurcation_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

result = struct('L6',l6,'Inhibition',inhibition,'Summary',summary);
l6ns_plot_bifurcations(result,outputDir);
save(fullfile(outputDir,'l6_inhibition_bifurcation_result.mat'),'result','-v7.3');
end

function out = local_existing_l6_result(cfg,pathway,l6Result)
t = l6Result.Continuation.Table;
if ~ismember('persistentMaxReal',t.Properties.VariableNames)
    t.persistentMaxReal = nan(height(t),1);
    for row = 1:height(t)
        gamma = 1-t.w(row);
        t.persistentMaxReal(row) = l6ns_max_real( ...
            l6ns_pathway_jacobian(pathway,gamma,1),cfg);
    end
end
t.freezeWeight = t.w;
t.dynamicGain = 1-t.w;
summary = l6Result.Summary;
out = struct();
out.Name = 'L6';
out.CriticalFreezeWeight = summary.wCritical(1);
out.CriticalGain = summary.alphaCritical(1);
out.Quadratic = summary.quadraticCoefficient(1);
out.Cubic = summary.cubicCoefficient(1);
out.Beta = summary.betaReal(1);
out.ZeroResidual = summary.compensatedZeroStateResidual(1);
out.StabilityExchange = logical(summary.stabilityExchange(1));
out.Classification = string(summary.classification(1));
out.Continuation = t;
end

function out = local_inhibition_result(cfg,pathway,context,fixed)
gammaGrid = primary_gamma_grid();
alpha = nan(size(gammaGrid));
for index = 1:numel(gammaGrid)
    alpha(index) = l6ns_max_real( ...
        l6ns_pathway_jacobian(pathway,1,gammaGrid(index)),cfg);
end
bracket = find((alpha(1:end-1)-1).*(alpha(2:end)-1)<=0,1,'last');
if isempty(bracket)
    error('No inhibition stability crossing was found in the search interval.');
end
low = gammaGrid(bracket);
high = gammaGrid(bracket+1);
fLow = alpha(bracket)-1;
for iteration = 1:34
    middle = (low+high)/2;
    value = l6ns_max_real(l6ns_pathway_jacobian(pathway,1,middle),cfg)-1;
    if fLow*value<=0
        high = middle;
    else
        low = middle;
        fLow = value;
    end
end
gammaCritical = (low+high)/2;
jCritical = l6ns_pathway_jacobian(pathway,1,gammaCritical);
modes = l6ns_eigenpairs(jCritical,10,'largestreal',cfg);
[~,criticalIndex] = min(abs(modes.Lambda-1));
r = real(modes.Right(:,criticalIndex));
r = r/norm(r);
l = real(modes.Left(:,criticalIndex));
l = l/conj(l'*r);
l = l/(l'*r);
beta = real(l'*pathway.JI*r);

phi = @(x,gamma)l6ns_phi(x,0,context,[1 1],gamma);
[normalForm,quadratic,cubic] = local_normal_form( ...
    phi,fixed,r,l,gammaCritical);
critical = struct('Parameter',gammaCritical,'Right',r,'Left',l, ...
    'Beta',beta,'Quadratic',quadratic,'Cubic',cubic);
continuation = local_continue_gain_branch(phi,fixed,critical,cfg);
continuation.persistentMaxReal = nan(height(continuation),1);
for row = 1:height(continuation)
    continuation.persistentMaxReal(row) = l6ns_max_real( ...
        l6ns_pathway_jacobian(pathway,1,continuation.dynamicGain(row)),cfg);
end

valid = continuation.converged & isfinite(continuation.maxRealLambda);
persistentStable = continuation.persistentMaxReal<1;
secondaryStable = continuation.maxRealLambda<1;
exchange = any(valid & persistentStable & ~secondaryStable) && ...
    any(valid & ~persistentStable & secondaryStable);
quadraticNonzero = abs(quadratic)>1e-5;
cubicNonzero = abs(cubic)>1e-8;
continuedBoth = any(valid & continuation.amplitude<0) && ...
    any(valid & continuation.amplitude>0);
if quadraticNonzero && exchange
    classification = "transcritical_supported";
elseif ~quadraticNonzero && cubicNonzero && continuedBoth && exchange
    classification = "pitchfork_like_supported";
elseif ~quadraticNonzero && cubicNonzero && continuedBoth
    classification = "pitchfork_or_symmetry_driven_candidate";
else
    classification = "inconclusive";
end

out = struct();
out.Name = 'inhibition';
out.CriticalFreezeWeight = 1-gammaCritical;
out.CriticalGain = gammaCritical;
out.Quadratic = quadratic;
out.Cubic = cubic;
out.Beta = beta;
out.ZeroResidual = norm(phi(zeros(size(fixed)),gammaCritical));
out.StabilityExchange = exchange;
out.Classification = classification;
out.Continuation = continuation;
out.NormalForm = normalForm;
end

function grid = primary_gamma_grid()
grid = linspace(0.72,1.02,31);
end

function [normalForm,quadraticEstimate,cubicEstimate] = local_normal_form( ...
        phi,fixed,r,l,gammaCritical)
stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
epsilon = stateScale*[1e-2 3e-3 1e-3 3e-4 1e-4];
quadratic = nan(size(epsilon));
cubic = nan(size(epsilon));
for index = 1:numel(epsilon)
    h = epsilon(index);
    center = phi(fixed,gammaCritical);
    plus = phi(fixed+h*r,gammaCritical);
    minus = phi(fixed-h*r,gammaCritical);
    quadratic(index) = 0.5*real(l'*(plus-2*center+minus))/(h^2);
    plus2 = phi(fixed+2*h*r,gammaCritical);
    minus2 = phi(fixed-2*h*r,gammaCritical);
    cubic(index) = real(l'*(plus2-2*plus+2*minus-minus2))/(12*h^3);
end
normalForm = table(epsilon(:),quadratic(:),cubic(:), ...
    'VariableNames',{'epsilon','quadraticCoefficient','cubicCoefficient'});
quadraticEstimate = local_plateau(quadratic);
cubicEstimate = local_plateau(cubic);
end

function estimate = local_plateau(values)
values = values(:);
score = nan(numel(values)-2,1);
medians = nan(size(score));
for index = 1:numel(score)
    window = values(index:index+2);
    medians(index) = median(window);
    score(index) = std(window)/max(abs(medians(index)),1e-12);
end
[~,best] = min(score);
estimate = medians(best);
end

function tableOut = local_continue_gain_branch(phi,fixed,critical,cfg)
r = critical.Right;
l = critical.Left;
scale = max(1,norm(fixed)/sqrt(numel(fixed)));
amplitudes = scale*[-0.025 -0.015 -0.008 -0.004 0.004 0.008 0.015 0.025];
rows = cell(numel(amplitudes),1);
for ai = 1:numel(amplitudes)
    target = amplitudes(ai);
    parameter = critical.Parameter-real((critical.Quadratic*target+ ...
        critical.Cubic*target^2)/critical.Beta);
    state = fixed+target*r;
    converged = false;
    residualNorm = inf;
    newtonIteration = 0;
    tolerance = 5e-10*max(1,norm(fixed));
    for newtonIteration = 1:18
        g = phi(state,parameter)-state;
        residual = [g;l'*(state-fixed)-target];
        residualNorm = norm(residual);
        if residualNorm<tolerance
            converged = true;
            break
        end
        stateStep = 5e-7*max(1,norm(state)/sqrt(numel(state)));
        parameterStep = 5e-7;
        phiBase = phi(state,parameter);
        dGdp = (phi(state,parameter+parameterStep)- ...
            phi(state,parameter-parameterStep))/(2*parameterStep);
        operator = @(delta)local_augmented_action(delta,state,parameter,phi, ...
            phiBase,stateStep,dGdp,l);
        [delta,~] = gmres(operator,-residual,[],1e-6,60);
        if any(~isfinite(delta)); break; end
        accepted = false;
        lineScale = 1;
        for lineIteration = 1:8
            candidateState = state+lineScale*delta(1:end-1);
            candidateParameter = parameter+lineScale*delta(end);
            candidateResidual = [phi(candidateState,candidateParameter)-candidateState; ...
                l'*(candidateState-fixed)-target];
            if norm(candidateResidual)<residualNorm
                state = candidateState;
                parameter = candidateParameter;
                accepted = true;
                break
            end
            lineScale = lineScale/2;
        end
        if ~accepted; break; end
    end
    maxReal = NaN;
    if converged
        nonlinearModes = l6ns_numeric_leading_modes( ...
            @(x)phi(x,parameter),state,4,cfg);
        maxReal = max(real(nonlinearModes.Lambda));
    end
    rows{ai} = {target,parameter,1-parameter,norm(state-fixed), ...
        residualNorm,residualNorm/max(1,norm(fixed)),converged, ...
        newtonIteration,maxReal};
    fprintf('Inhibition continuation a=%+.5g gamma=%+.9g residual %.3e converged=%d.\n', ...
        target,parameter,residualNorm,converged);
end
tableOut = cell2table(vertcat(rows{:}),'VariableNames', ...
    {'amplitude','dynamicGain','freezeWeight','branchDistance','residualNorm', ...
    'relativeResidual','converged','newtonIterations','maxRealLambda'});
end

function output = local_augmented_action(delta,state,parameter,phi,phiBase,stateStep,dGdp,l)
direction = delta(1:end-1);
parameterDirection = delta(end);
directionNorm = norm(direction);
if directionNorm==0
    dGdf = zeros(size(direction));
else
    h = stateStep/directionNorm;
    dGdf = (phi(state+h*direction,parameter)-phiBase)/h-direction;
end
output = [dGdf+dGdp*parameterDirection;l'*direction];
end

function tableOut = local_stability_curve(cfg,pathway,freezeWeight,pathwayName)
alpha = nan(size(freezeWeight));
for index = 1:numel(freezeWeight)
    gain = 1-freezeWeight(index);
    if strcmp(pathwayName,'L6')
        matrix = l6ns_pathway_jacobian(pathway,gain,1);
    else
        matrix = l6ns_pathway_jacobian(pathway,1,gain);
    end
    alpha(index) = l6ns_max_real(matrix,cfg);
end
tableOut = table(freezeWeight(:),1-freezeWeight(:),alpha(:),alpha(:)<1, ...
    'VariableNames',{'freezeWeight','dynamicGain','maxRealLambda','stable'});
end
