function diagnostic = l6ns_atlas_point_metrics(setup,pathwayName,state,w,dwds,options)
% Compute equilibrium, spectrum, singularity, population, and clamp metrics.

if nargin<5 || isempty(dwds); dwds=NaN; end
if nargin<6; options=struct(); end
modeCount = local_option(options,'ModeCount',40);
maximumModeCount = local_option(options,'MaximumModeCount',160);
boundaryTolerance = local_option(options,'BoundaryTolerance',1e-8);

state = state(:);
fixed = setup.Context.FixedPoint(:);
phi = local_phi(setup,pathwayName,w);
residual = phi(state)-state;
J = l6ns_state_jacobian(setup,pathwayName,state,w);

count = modeCount;
while true
    modes = l6ns_eigenpairs(J,count,'largestreal',setup.Config);
    if count>=maximumModeCount || real(modes.Lambda(end))<1-boundaryTolerance
        break
    end
    count = min(maximumModeCount,2*count);
end
unstableDimension = sum(real(modes.Lambda)>1+boundaryTolerance);
nearBoundary = abs(real(modes.Lambda)-1)<=boundaryTolerance;
criticalMultiplicity = sum(nearBoundary);

A = J-speye(size(J,1));
singularValues = [NaN;NaN];
try
    singularValues = sort(svds(A,2,'smallest'));
catch
    try
        singularValues = sort(svds(A,2,0));
    catch
    end
end

n = numel(state)/3;
population = reshape(state,n,3);
fixedPopulation = reshape(fixed,n,3);
populationMean = mean(population,1);
populationMinimum = min(population,[],1);
populationMaximum = max(population,[],1);
populationRmsDistance = sqrt(mean((population-fixedPopulation).^2,1));
negativeFraction = mean(population<0,1);
fixedMean = mean(fixedPopulation,1);
meanRateRatio = populationMean./max(fixedMean,eps);
maximumRateRatio = max(state)/max(max(fixed),eps);
extremeHighRate = maximumRateRatio>=2 || any(meanRateRatio>=2);
[lowerMargin,upperMargin,activeFraction] = local_l6_clamp_metrics(state,setup.Context);

[~,boundaryOrder] = sort(abs(real(modes.Lambda)-1),'ascend');
boundaryOrder = boundaryOrder(1:min(modeCount,numel(boundaryOrder)));
modeTable = table((1:numel(boundaryOrder))',modes.Lambda(boundaryOrder), ...
    modes.ConditionNumber(boundaryOrder),modes.RightResidual(boundaryOrder), ...
    modes.LeftResidual(boundaryOrder), ...
    'VariableNames',{'boundaryRank','lambda','conditionNumber', ...
    'rightResidual','leftResidual'});

summary = table(string(pathwayName),w,norm(residual), ...
    norm(residual)/max(1,norm(state)),max(real(modes.Lambda)), ...
    unstableDimension,criticalMultiplicity,singularValues(1),singularValues(2),dwds, ...
    norm(state-fixed),min(state),max(state), ...
    populationMean(1),populationMean(2),populationMean(3), ...
    populationMinimum(1),populationMinimum(2),populationMinimum(3), ...
    populationMaximum(1),populationMaximum(2),populationMaximum(3), ...
    populationRmsDistance(1),populationRmsDistance(2),populationRmsDistance(3), ...
    negativeFraction(1),negativeFraction(2),negativeFraction(3), ...
    meanRateRatio(1),meanRateRatio(2),meanRateRatio(3),maximumRateRatio, ...
    extremeHighRate,lowerMargin,upperMargin,activeFraction, ...
    'VariableNames',{'pathway','freezeWeight','residualNorm','relativeResidual', ...
    'maxRealLambda','unstableDimension','criticalMultiplicity','sigmaMinA','sigma2A', ...
    'dwds','distanceFromPersistent','minimumState','maximumState', ...
    'meanS','meanC','meanI','minimumS','minimumC','minimumI', ...
    'maximumS','maximumC','maximumI','rmsDistanceS','rmsDistanceC','rmsDistanceI', ...
    'negativeFractionS','negativeFractionC','negativeFractionI', ...
    'meanRateRatioS','meanRateRatioC','meanRateRatioI','maximumRateRatio', ...
    'extremeHighRate', ...
    'lowerClampMargin','upperClampMargin','clampActiveFraction'});

diagnostic = struct('Summary',summary,'Modes',modeTable, ...
    'RightVectors',modes.Right(:,boundaryOrder), ...
    'LeftVectors',modes.Left(:,boundaryOrder),'Jacobian',J);
end

function phi = local_phi(setup,pathwayName,w)
context=setup.Context;
if strcmpi(pathwayName,'L6')
    phi=@(x)l6ns_phi(x,w,context,[1 1],1);
elseif any(strcmpi(pathwayName,{'I','Inhibition'}))
    phi=@(x)l6ns_phi(x,0,context,[1 1],1-w);
else
    error('Unknown pathway %s.',pathwayName);
end
end

function [lowerMargin,upperMargin,activeFraction] = local_l6_clamp_metrics(state,context)
[~,raw] = l6ns_l6_dynamic(state,context,[1 1]);
lowerMargin = min(raw-1);
upperMargin = min(40-raw);
activeFraction = mean(raw<=1 | raw>=40);
end

function value=local_option(options,name,defaultValue)
if isfield(options,name); value=options.(name); else; value=defaultValue; end
end
