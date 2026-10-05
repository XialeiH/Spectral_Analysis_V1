function setup = l6ns_global_bifurcation_setup(cfg,context,pathway)
% Build the reusable context and critical data for wide branch continuation.

outputDir = fullfile(cfg.OutputRoot,'global_bifurcation_extension');
if ~exist(outputDir,'dir'); mkdir(outputDir); end

fixed = context.FixedPoint(:);
rawContext = context;
rawContext.FixedPointCorrection = [];
context.FixedPointCorrection = fixed-l6ns_phi(fixed,0,rawContext,[1 1],1);

l6Critical = local_critical_data(cfg,context,pathway,'L6');
iCritical = local_critical_data(cfg,context,pathway,'I');
symmetry = local_symmetry_audit(context,l6Critical,iCritical,cfg.RandomSeed);

writetable(symmetry,fullfile(outputDir,'spatial_symmetry_audit.tsv'), ...
    'FileType','text','Delimiter','\t');

setup = struct();
setup.Config = cfg;
setup.Context = context;
setup.Pathway = pathway;
setup.L6 = l6Critical;
setup.Inhibition = iCritical;
setup.SymmetryAudit = symmetry;
setup.OutputDir = outputDir;
save(fullfile(outputDir,'global_bifurcation_setup.mat'),'setup','-v7.3');
end

function critical = local_critical_data(cfg,context,pathway,pathwayName)
if strcmp(pathwayName,'L6')
    matrixAt = @(w)l6ns_pathway_jacobian(pathway,1-w,1);
    phiAt = @(x,w)l6ns_phi(x,w,context,[1 1],1);
    bracket = [-0.30 -0.05];
    derivativeMatrix = -pathway.J6;
else
    matrixAt = @(w)l6ns_pathway_jacobian(pathway,1,1-w);
    phiAt = @(x,w)l6ns_phi(x,0,context,[1 1],1-w);
    bracket = [0.05 0.25];
    derivativeMatrix = -pathway.JI;
end

low = bracket(1);
high = bracket(2);
fLow = l6ns_max_real(matrixAt(low),cfg)-1;
fHigh = l6ns_max_real(matrixAt(high),cfg)-1;
if fLow*fHigh>0
    error('%s critical point is not bracketed by [%.3g, %.3g].', ...
        pathwayName,low,high);
end
for iteration = 1:42
    middle = (low+high)/2;
    fMiddle = l6ns_max_real(matrixAt(middle),cfg)-1;
    if fLow*fMiddle<=0
        high = middle;
        fHigh = fMiddle; %#ok<NASGU>
    else
        low = middle;
        fLow = fMiddle;
    end
end
wc = (low+high)/2;
jCritical = matrixAt(wc);
modes = l6ns_eigenpairs(jCritical,20,'largestreal',cfg);
[~,criticalIndex] = min(abs(modes.Lambda-1));
r = real(modes.Right(:,criticalIndex));
r = r/norm(r);
l = real(modes.Left(:,criticalIndex));
l = l/(l'*r);

other = modes.Lambda;
other(criticalIndex) = [];
beta = real(l'*derivativeMatrix*r);
[normalForm,quadratic,cubic] = local_normal_form(phiAt,context.FixedPoint,r,l,wc);
writetable(normalForm,fullfile(cfg.OutputRoot,'global_bifurcation_extension', ...
    sprintf('%s_normal_form.tsv',lower(pathwayName))), ...
    'FileType','text','Delimiter','\t');

critical = struct();
critical.Name = pathwayName;
critical.W = wc;
critical.Right = r;
critical.Left = l;
critical.Beta = beta;
critical.Quadratic = quadratic;
critical.Cubic = cubic;
critical.Lambda = modes.Lambda(criticalIndex);
critical.SpectralSeparation = min(abs(other-modes.Lambda(criticalIndex)));
critical.SecondLeadingReal = max(real(other));
critical.NormalForm = normalForm;
end

function [tableOut,quadraticEstimate,cubicEstimate] = local_normal_form( ...
        phi,fixed,r,l,wc)
stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
epsilon = stateScale*[1e-2 3e-3 1e-3 3e-4 1e-4];
quadratic = nan(size(epsilon));
cubic = nan(size(epsilon));
center = phi(fixed,wc);
for index = 1:numel(epsilon)
    h = epsilon(index);
    plus = phi(fixed+h*r,wc);
    minus = phi(fixed-h*r,wc);
    quadratic(index) = 0.5*real(l'*(plus-2*center+minus))/(h^2);
    plus2 = phi(fixed+2*h*r,wc);
    minus2 = phi(fixed-2*h*r,wc);
    cubic(index) = real(l'*(plus2-2*plus+2*minus-minus2))/(12*h^3);
end
tableOut = table(epsilon(:),quadratic(:),cubic(:), ...
    'VariableNames',{'epsilon','quadraticCoefficient','cubicCoefficient'});
quadraticEstimate = local_plateau(quadratic);
cubicEstimate = local_plateau(cubic);
end

function estimate = local_plateau(values)
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

function tableOut = local_symmetry_audit(context,l6Critical,iCritical,seed)
mapSize = context.MapSize;
if any(mod(mapSize,4)~=0)
    error('Map size [%d %d] is not divisible by four.',mapSize(1),mapSize(2));
end
fixed = context.FixedPoint(:);
stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
rng(seed);
probe = fixed+1e-3*stateScale*randn(size(fixed));
rows = {};
for rowShift = 0:mapSize(1)/4:mapSize(1)-mapSize(1)/4
    for columnShift = 0:mapSize(2)/4:mapSize(2)-mapSize(2)/4
        shiftedFixed = local_shift_state(fixed,mapSize,rowShift,columnShift);
        fixedError = norm(shiftedFixed-fixed)/max(norm(fixed),eps);
        for pathwayIndex = 1:2
            if pathwayIndex==1
                critical = l6Critical;
                phi = @(x)l6ns_phi(x,critical.W,context,[1 1],1);
            else
                critical = iCritical;
                phi = @(x)l6ns_phi(x,0,context,[1 1],1-critical.W);
            end
            shiftedRight = local_shift_state(critical.Right,mapSize,rowShift,columnShift);
            plusError = norm(shiftedRight-critical.Right)/norm(critical.Right);
            minusError = norm(shiftedRight+critical.Right)/norm(critical.Right);
            shiftedProbe = local_shift_state(probe,mapSize,rowShift,columnShift);
            equivarianceError = norm(phi(shiftedProbe)- ...
                local_shift_state(phi(probe),mapSize,rowShift,columnShift)) / ...
                max(norm(phi(probe)),eps);
            rows(end+1,:) = {string(critical.Name),rowShift,columnShift, ...
                fixedError,plusError,minusError,equivarianceError}; %#ok<AGROW>
        end
    end
end
tableOut = cell2table(rows,'VariableNames',{'pathway','rowShift','columnShift', ...
    'fixedPointInvarianceError','criticalModeEvenError','criticalModeOddError', ...
    'mapEquivarianceError'});
end

function shifted = local_shift_state(state,mapSize,rowShift,columnShift)
n = prod(mapSize);
shifted = zeros(size(state));
for population = 1:3
    indices = (population-1)*n+(1:n);
    map = reshape(state(indices),mapSize);
    shifted(indices) = reshape(circshift(map,[rowShift columnShift]),[],1);
end
end
