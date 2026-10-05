function branch = l6ns_pathway_branch(pathway,gammaGrid,axisName,cfg)
% Track one biorthogonal branch while separately retaining spectral abscissa.

nPoint = numel(gammaGrid);
n = pathway.PopulationSize;
modePools = cell(nPoint,1);
for gi = 1:nPoint
    [gamma6,gammaI] = local_gains(gammaGrid(gi),axisName);
    jacobian = l6ns_pathway_jacobian(pathway,gamma6,gammaI);
    modePools{gi} = l6ns_eigenpairs( ...
        jacobian,cfg.BranchPoolSize,'largestreal',cfg);
end

[~,anchor] = min(abs(gammaGrid-1));
selected = nan(nPoint,1);
[~,selected(anchor)] = max(real(modePools{anchor}.Lambda));
for gi = anchor+1:nPoint
    selected(gi) = local_match(modePools{gi-1},selected(gi-1),modePools{gi});
end
for gi = anchor-1:-1:1
    selected(gi) = local_match(modePools{gi+1},selected(gi+1),modePools{gi});
end

rows = cell(nPoint,1);
right = cell(nPoint,1);
left = cell(nPoint,1);
for gi = 1:nPoint
    gamma = gammaGrid(gi);
    [gamma6,gammaI] = local_gains(gamma,axisName);
    jacobian = l6ns_pathway_jacobian(pathway,gamma6,gammaI);
    modes = modePools{gi};
    index = selected(gi);
    lambda = modes.Lambda(index);
    spectralAbscissa = max(real(modes.Lambda));
    r = modes.Right(:,index);
    l = modes.Left(:,index);
    right{gi} = r;
    left{gi} = l;

    if gi == anchor
        continuity = 1;
    elseif gi < anchor
        continuity = local_pair_overlap(l,r,modePools{gi+1},selected(gi+1));
    else
        continuity = local_pair_overlap(l,r,modePools{gi-1},selected(gi-1));
    end

    generator = (jacobian-speye(size(jacobian,1))) / cfg.TauMs;
    hermitianPart = (generator+generator')/2;
    numericalAbscissa = real(eigs(hermitianPart,1,'largestreal'));
    stabilityMargin = 1-spectralAbscissa;
    recoveryMs = Inf;
    if stabilityMargin > 0
        recoveryMs = cfg.TauMs/stabilityMargin;
    end
    frequencyHz = abs(imag(lambda))*1000/(2*pi*cfg.TauMs);
    s6 = l'*pathway.J6*r;
    sI = l'*pathway.JI*r;
    target6 = local_target_terms(pathway.J6,l,r,n);
    targetI = local_target_terms(pathway.JI,l,r,n);
    closure = l'*jacobian*r-lambda;
    rows{gi} = {gamma,gamma6,gammaI,spectralAbscissa,real(lambda),imag(lambda), ...
        (spectralAbscissa-1)/cfg.TauMs,stabilityMargin,recoveryMs,frequencyHz, ...
        modes.ConditionNumber(index),1/modes.ConditionNumber(index), ...
        numericalAbscissa,real(s6),imag(s6),real(sI),imag(sI), ...
        real(target6(1)),real(target6(2)),real(target6(3)), ...
        real(targetI(1)),real(targetI(2)),real(targetI(3)), ...
        abs(closure),modes.RightResidual(index),modes.LeftResidual(index),continuity};
end

tableOut = cell2table(vertcat(rows{:}),'VariableNames',{ ...
    'gamma','gamma6','gammaI','spectralAbscissaJ','lambdaReal','lambdaImag', ...
    'alphaPhysicalPerMs','stabilityMarginJ','recoveryTimeMs','frequencyHz', ...
    'conditionNumber','leftRightAlignment','numericalAbscissaPerMs', ...
    'sensitivity6Real','sensitivity6Imag','sensitivityIReal','sensitivityIImag', ...
    'L6toSReal','L6toCReal','L6toIReal','ItoSReal','ItoCReal','ItoIReal', ...
    'attributionClosure','rightResidual','leftResidual','adjacentModeOverlap'});

if nPoint >= 3
    finiteDifference = gradient(tableOut.lambdaReal,tableOut.gamma);
else
    finiteDifference = nan(nPoint,1);
end
if strcmp(axisName,'L6')
    analytic = tableOut.sensitivity6Real;
else
    analytic = tableOut.sensitivityIReal;
end
tableOut.finiteDifferenceDerivative = finiteDifference;
tableOut.derivativeRelativeError = abs(finiteDifference-analytic) ./ ...
    max(abs(analytic),1e-12);

branch = struct('Axis',axisName,'Table',tableOut,'Right',{right}, ...
    'Left',{left},'SelectedIndex',selected);
end

function [gamma6,gammaI] = local_gains(gamma,axisName)
if strcmp(axisName,'L6')
    gamma6 = gamma;
    gammaI = 1;
elseif strcmp(axisName,'I')
    gamma6 = 1;
    gammaI = gamma;
else
    error('axisName must be L6 or I.');
end
end

function index = local_match(previousModes,previousIndex,currentModes)
previousLambda = previousModes.Lambda(previousIndex);
previousRight = previousModes.Right(:,previousIndex);
previousLeft = previousModes.Left(:,previousIndex);
score = nan(numel(currentModes.Lambda),1);
scale = max(0.03,0.12*max(1,abs(previousLambda)));
for candidate = 1:numel(currentModes.Lambda)
    forward = abs(previousLeft'*currentModes.Right(:,candidate)) / ...
        max(norm(previousLeft)*norm(currentModes.Right(:,candidate)),eps);
    backward = abs(currentModes.Left(:,candidate)'*previousRight) / ...
        max(norm(currentModes.Left(:,candidate))*norm(previousRight),eps);
    proximity = exp(-abs(currentModes.Lambda(candidate)-previousLambda)/scale);
    score(candidate) = 0.4*forward+0.4*backward+0.2*proximity;
end
[~,index] = max(score);
end

function overlap = local_pair_overlap(l,r,otherModes,otherIndex)
otherL = otherModes.Left(:,otherIndex);
otherR = otherModes.Right(:,otherIndex);
forward = abs(l'*otherR)/max(norm(l)*norm(otherR),eps);
backward = abs(otherL'*r)/max(norm(otherL)*norm(r),eps);
overlap = sqrt(forward*backward);
end

function terms = local_target_terms(matrix,l,r,n)
terms = zeros(3,1);
for population = 1:3
    rows = (population-1)*n+(1:n);
    terms(population) = l(rows)'*(matrix(rows,:)*r);
end
end
