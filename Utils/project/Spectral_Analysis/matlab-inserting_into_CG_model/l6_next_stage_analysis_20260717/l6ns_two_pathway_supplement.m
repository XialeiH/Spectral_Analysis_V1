function supplement = l6ns_two_pathway_supplement(cfg,data,context)
% Complete the joint reduced, synthetic-readout, and boundary deliverables.

inputFile = fullfile(cfg.OutputRoot,'two_pathway_plan','two_pathway_result.mat');
loaded = load(inputFile,'result');
primary = loaded.result;
outputDir = fullfile(cfg.OutputRoot,'two_pathway_supplement');
if ~exist(outputDir,'dir'); mkdir(outputDir); end

reduced = local_reduced_model(cfg,primary,outputDir);
synthetic = local_synthetic_detectability(primary,outputDir);
boundary = local_inhibition_boundary(cfg,primary,context,data,outputDir);
supplement = struct('Reduced',reduced,'Synthetic',synthetic,'Boundary',boundary);
save(fullfile(outputDir,'two_pathway_supplement_result.mat'),'supplement','-v7.3');
end

function reduced = local_reduced_model(cfg,primary,outputDir)
pathway = primary.Pathway;
weights = primary.TransientAudit.Physical.WeightVector(:);
inverseWeights = 1./weights;

candidate = [];
candidate = local_add_branch(candidate,primary.L6Branch,weights,inverseWeights);
candidate = local_add_branch(candidate,primary.InhibitionBranch,weights,inverseWeights);
for pairIndex = 1:numel(primary.Matched.Pairs)
    pair = primary.Matched.Pairs{pairIndex};
    candidate = [candidate, ...
        weights.*pair.Modes6.Right(:,1:4), ...
        inverseWeights.*pair.Modes6.Left(:,1:4), ...
        weights.*pair.ModesI.Right(:,1:4), ...
        inverseWeights.*pair.ModesI.Left(:,1:4), ...
        weights.*pair.Physical6.OptimalInput, ...
        weights.*pair.Physical6.PeakOutput, ...
        weights.*pair.PhysicalI.OptimalInput, ...
        weights.*pair.PhysicalI.PeakOutput]; %#ok<AGROW>
end

gridTable = primary.Interaction.Table;
gridModes = cell(height(gridTable),1);
for row = 1:height(gridTable)
    fullJ = l6ns_pathway_jacobian(pathway, ...
        gridTable.gamma6(row),gridTable.gammaI(row));
    gridModes{row} = l6ns_eigenpairs(fullJ,4,'largestreal',cfg);
    candidate = [candidate, ...
        weights.*gridModes{row}.Right(:,1:4), ...
        inverseWeights.*gridModes{row}.Left(:,1:4)]; %#ok<AGROW>
end

scaledBaseline = spdiags(weights,0,numel(weights),numel(weights)) * ...
    pathway.JBaseline * spdiags(inverseWeights,0,numel(weights),numel(weights));
controlKrylov = weights.*primary.Operators.B;
observationKrylov = inverseWeights.*primary.Operators.C';
for order = 0:3
    candidate = [candidate,controlKrylov,observationKrylov]; %#ok<AGROW>
    controlKrylov = scaledBaseline*controlKrylov;
    observationKrylov = scaledBaseline'*observationKrylov;
end
candidate = full(candidate);
candidateNorm = sqrt(sum(abs(candidate).^2,1));
candidate = candidate(:,candidateNorm>eps);
candidate = candidate ./ sqrt(sum(abs(candidate).^2,1));
[leftCandidate,singularCandidate,~] = svd(candidate,'econ');
candidateSpectrum = diag(singularCandidate);
retain = min(128,sum(candidateSpectrum>1e-10*candidateSpectrum(1)));
q = leftCandidate(:,1:retain);

jRest = q'*(weights.*(pathway.JRest*(inverseWeights.*q)));
j6 = q'*(weights.*(pathway.J6*(inverseWeights.*q)));
jI = q'*(weights.*(pathway.JI*(inverseWeights.*q)));
bReduced = q'*(weights.*primary.Operators.B);
cReduced = primary.Operators.C*(inverseWeights.*q);

rows = cell(height(gridTable),1);
for row = 1:height(gridTable)
    gamma6 = gridTable.gamma6(row);
    gammaI = gridTable.gammaI(row);
    reducedJ = jRest+gamma6*j6+gammaI*jI;
    [vr,dr] = eig(full(reducedJ),'vector');
    [~,critical] = max(real(dr));
    reducedLambda = dr(critical);
    fullModes = gridModes{row};
    [~,fullCritical] = max(real(fullModes.Lambda));
    rightFull = weights.*fullModes.Right(:,fullCritical);
    rightReduced = q*vr(:,critical);
    rightCosine = abs(rightFull'*rightReduced) / ...
        max(norm(rightFull)*norm(rightReduced),eps);

    [vl,dl] = eig(full(reducedJ'),'vector');
    [~,leftIndex] = min(abs(conj(dl)-reducedLambda));
    leftFull = inverseWeights.*fullModes.Left(:,fullCritical);
    leftReduced = q*vl(:,leftIndex);
    leftCosine = abs(leftFull'*leftReduced) / ...
        max(norm(leftFull)*norm(leftReduced),eps);

    [physicalGain,constrainedGain] = local_reduced_gains( ...
        reducedJ,bReduced,cReduced,cfg.GridTimeGridMs,cfg.TauMs);
    rows{row} = {gamma6,gammaI,gridTable.lambdaMaxReal(row),real(reducedLambda), ...
        abs(real(reducedLambda)-gridTable.lambdaMaxReal(row)),rightCosine,leftCosine, ...
        gridTable.physicalPeakGain(row),physicalGain, ...
        abs(physicalGain-gridTable.physicalPeakGain(row))/ ...
            max(gridTable.physicalPeakGain(row),eps), ...
        gridTable.constrainedPeakGain(row),constrainedGain, ...
        abs(constrainedGain-gridTable.constrainedPeakGain(row))/ ...
            max(gridTable.constrainedPeakGain(row),eps)};
end
validation = cell2table(vertcat(rows{:}),'VariableNames', ...
    {'gamma6','gammaI','fullLambdaReal','reducedLambdaReal','lambdaAbsoluteError', ...
    'rightModeCosine','leftModeCosine','fullPhysicalGain','reducedPhysicalGain', ...
    'physicalGainRelativeError','fullConstrainedGain','reducedConstrainedGain', ...
    'constrainedGainRelativeError'});
writetable(validation,fullfile(outputDir,'joint_reduced_validation.tsv'), ...
    'FileType','text','Delimiter','\t');

baselineJ = jRest+j6+jI;
[r,l,lambda] = local_reduced_critical_pair(baselineJ);
s6 = real(l'*j6*r);
sI = real(l'*jI*r);
reducedSlope = -s6/sI;
fullSlope = primary.Interaction.Compensation.predictedDGammaIDGamma6;
summary = table(retain,max(validation.lambdaAbsoluteError), ...
    median(validation.lambdaAbsoluteError),min(validation.rightModeCosine), ...
    min(validation.leftModeCosine),max(validation.physicalGainRelativeError), ...
    median(validation.physicalGainRelativeError), ...
    max(validation.constrainedGainRelativeError), ...
    median(validation.constrainedGainRelativeError),real(lambda), ...
    s6,sI,fullSlope,reducedSlope,abs(reducedSlope-fullSlope), ...
    'VariableNames',{'dimension','maximumLambdaError','medianLambdaError', ...
    'minimumRightModeCosine','minimumLeftModeCosine', ...
    'maximumPhysicalGainRelativeError','medianPhysicalGainRelativeError', ...
    'maximumConstrainedGainRelativeError','medianConstrainedGainRelativeError', ...
    'baselineReducedLambdaReal','reducedS6','reducedSI','fullCompensationSlope', ...
    'reducedCompensationSlope','compensationSlopeAbsoluteError'});
writetable(summary,fullfile(outputDir,'joint_reduced_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[80 80 1200 420]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; scatter(validation.fullLambdaReal,validation.reducedLambdaReal,40,'filled');
hold on; limits=xlim; plot(limits,limits,'k--'); grid on; axis square;
xlabel('full max Re \lambda'); ylabel('reduced max Re \lambda');
nexttile; scatter(validation.fullPhysicalGain,validation.reducedPhysicalGain,40,'filled');
hold on; limits=xlim; plot(limits,limits,'k--'); grid on; axis square;
xlabel('full physical gain'); ylabel('reduced physical gain');
nexttile; scatter(validation.fullConstrainedGain,validation.reducedConstrainedGain,40,'filled');
hold on; limits=xlim; plot(limits,limits,'k--'); grid on; axis square;
xlabel('full constrained gain'); ylabel('reduced constrained gain');
l6ns_save_figure(fig,outputDir,'figure6_joint_reduced_validation'); close(fig);
reduced = struct('BasisScaled',q,'CandidateSingularValues',candidateSpectrum, ...
    'JRest',jRest,'J6',j6,'JI',jI,'Validation',validation,'Summary',summary);
end

function candidate = local_add_branch(candidate,branch,weights,inverseWeights)
indices = unique([1,round(numel(branch.Right)/2),numel(branch.Right)]);
for index = indices
    candidate = [candidate,weights.*branch.Right{index}, ...
        inverseWeights.*branch.Left{index}]; %#ok<AGROW>
end
end

function [physicalGain,constrainedGain] = local_reduced_gains(j,b,c,timeGrid,tau)
generator = (j-eye(size(j)))/tau;
physicalGain = 0;
constrainedGain = 0;
for time = timeGrid(:)'
    propagator = expm(time*generator);
    physicalGain = max(physicalGain,svds(propagator,1));
    if time>0
        constrainedGain = max(constrainedGain,svds(c*propagator*b,1));
    end
end
end

function [r,l,lambda] = local_reduced_critical_pair(j)
[right,values] = eig(full(j),'vector');
[~,index] = max(real(values));
lambda = values(index);
r = right(:,index);
[left,leftValues] = eig(full(j'),'vector');
[~,leftIndex] = min(abs(conj(leftValues)-lambda));
l = left(:,leftIndex);
l = l/conj(l'*r);
end

function synthetic = local_synthetic_detectability(primary,outputDir)
pair = primary.Matched.Pairs{end};
time = pair.Constrained6.TimeGridMs(:);
timeBin = (0:5:max(time))';
k6 = pair.Constrained6.Kernels;
kI = pair.ConstrainedI.Kernels;
observations = primary.Operators.ObservationNames;
controls = primary.Operators.ControlNames;
metricRows = {};
detectRows = {};
curveRows = {};
trialCount = 50;
noiseFraction = 0.05;
for observation = 1:numel(observations)
    for control = 1:numel(controls)
        curve6 = interp1(time,squeeze(k6(observation,control,:)),timeBin,'pchip');
        curveI = interp1(time,squeeze(kI(observation,control,:)),timeBin,'pchip');
        metrics6 = local_kernel_metrics(timeBin,curve6);
        metricsI = local_kernel_metrics(timeBin,curveI);
        metricRows(end+1,:) = [{observations{observation},controls{control},'L6'},metrics6]; %#ok<AGROW>
        metricRows(end+1,:) = [{observations{observation},controls{control},'I'},metricsI]; %#ok<AGROW>
        noiseSigma = noiseFraction*max(abs([curve6;curveI]));
        meanSe = noiseSigma/sqrt(trialCount);
        differenceSe = sqrt(2)*meanSe;
        aggregateZ = norm(curve6-curveI)/max(differenceSe,eps);
        maxBinZ = max(abs(curve6-curveI))/max(differenceSe,eps);
        detectRows(end+1,:) = {observations{observation},controls{control}, ...
            trialCount,noiseFraction,noiseSigma,aggregateZ,maxBinZ}; %#ok<AGROW>
        for ti = 1:numel(timeBin)
            curveRows(end+1,:) = {observations{observation},controls{control}, ...
                timeBin(ti),curve6(ti),curveI(ti),1.96*meanSe}; %#ok<AGROW>
        end
    end
end
metricTable = cell2table(metricRows,'VariableNames', ...
    {'observation','control','pathway','signedPeak','peakTimeMs','rebound', ...
    'tailTauMs','lateToEarlyArea','phase10HzRad'});
detectability = cell2table(detectRows,'VariableNames', ...
    {'observation','control','trials','noiseFraction','noiseSigma', ...
    'aggregateDifferenceZ','maximumBinDifferenceZ'});
curves = cell2table(curveRows,'VariableNames', ...
    {'observation','control','timeMs','l6Mean','inhibitionMean','mean95HalfWidth'});
writetable(metricTable,fullfile(outputDir,'synthetic_kernel_metrics.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(detectability,fullfile(outputDir,'synthetic_detectability.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(curves,fullfile(outputDir,'synthetic_binned_curves.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[80 80 1250 720]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
controlIndex = [1 3 5];
for observation = 1:2
    for column = 1:3
        nexttile;
        one = curves(strcmp(curves.observation,observations{observation}) & ...
            strcmp(curves.control,controls{controlIndex(column)}),:);
        plot(one.timeMs,one.l6Mean,'-','LineWidth',1.5); hold on;
        plot(one.timeMs,one.inhibitionMean,'--','LineWidth',1.5);
        plot(one.timeMs,one.l6Mean+one.mean95HalfWidth,':','Color',[0.4 0.4 0.4]);
        plot(one.timeMs,one.l6Mean-one.mean95HalfWidth,':','Color',[0.4 0.4 0.4]);
        grid on; xlabel('time (ms)'); ylabel('binned response');
        title(sprintf('%s, %s',observations{observation}, ...
            controls{controlIndex(column)}),'Interpreter','none');
        if observation==1 && column==1
            legend({'L6 matched','I matched','95% mean CI'},'Location','best');
        end
    end
end
l6ns_save_figure(fig,outputDir,'figure6_synthetic_detectability'); close(fig);
synthetic = struct('Metrics',metricTable,'Detectability',detectability,'Curves',curves);
end

function metrics = local_kernel_metrics(time,curve)
[~,peakIndex] = max(abs(curve));
signedPeak = curve(peakIndex);
peakTime = time(peakIndex);
later = curve(peakIndex:end);
if signedPeak>=0; rebound = min(later); else; rebound = max(later); end
early = trapz(time(time<=50),abs(curve(time<=50)));
late = trapz(time(time>=50),abs(curve(time>=50)));
lateToEarly = late/max(early,eps);
tailMask = time>=max(50,peakTime) & abs(curve)>max(abs(curve))*1e-6;
tailTau = nan;
if nnz(tailMask)>=4
    coefficients = polyfit(time(tailMask),log(abs(curve(tailMask))),1);
    if coefficients(1)<0; tailTau=-1/coefficients(1); end
end
frequency = 10/1000;
fourier = trapz(time,curve.*exp(-1i*2*pi*frequency*time));
metrics = {signedPeak,peakTime,rebound,tailTau,lateToEarly,angle(fourier)};
end

function boundary = local_inhibition_boundary(cfg,primary,context,data,outputDir)
branch = primary.InhibitionBranch.Table;
crossing = find((branch.spectralAbscissaJ(1:end-1)-1).* ...
    (branch.spectralAbscissaJ(2:end)-1)<=0,1,'last');
if isempty(crossing)
    boundary = struct('Status','no_boundary_in_sweep');
    return
end
low = branch.gamma(crossing); high = branch.gamma(crossing+1);
fLow = branch.spectralAbscissaJ(crossing)-1;
for iteration = 1:28
    middle = (low+high)/2;
    value = l6ns_max_real(l6ns_pathway_jacobian(primary.Pathway,1,middle),cfg)-1;
    if fLow*value<=0; high=middle; else; low=middle; fLow=value; end
end
gammaCritical = (low+high)/2;
jCritical = l6ns_pathway_jacobian(primary.Pathway,1,gammaCritical);
modes = l6ns_eigenpairs(jCritical,8,'largestreal',cfg);
[~,criticalIndex] = max(real(modes.Lambda));
lambda = modes.Lambda(criticalIndex);
r = real(modes.Right(:,criticalIndex)); r=r/norm(r);
l = real(modes.Left(:,criticalIndex)); l=l/conj(l'*r);
other = modes.Lambda; other(criticalIndex)=[];
separation = min(abs(other-lambda));
beta = l'*primary.Pathway.JI*r;

fixed = context.FixedPoint(:);
raw = context; raw.FixedPointCorrection=[];
context.FixedPointCorrection = fixed-l6ns_phi(fixed,0,raw,[1 1],1);
phi = @(x)l6ns_phi(x,0,context,[1 1],gammaCritical);
stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
epsilon = stateScale*[1e-2 3e-3 1e-3 3e-4 1e-4];
quadratic = nan(size(epsilon)); cubic = nan(size(epsilon));
phi0 = phi(fixed);
for index = 1:numel(epsilon)
    h = epsilon(index);
    plus = phi(fixed+h*r); minus = phi(fixed-h*r);
    quadratic(index) = 0.5*l'*(plus-2*phi0+minus)/(h^2);
    plus2 = phi(fixed+2*h*r); minus2 = phi(fixed-2*h*r);
    cubic(index) = (1/6)*l'*(plus2-2*plus+2*minus-minus2)/(2*h^3);
end
normalForm = table(epsilon(:),quadratic(:),cubic(:), ...
    'VariableNames',{'epsilon','quadraticCoefficient','cubicCoefficient'});
writetable(normalForm,fullfile(outputDir,'inhibition_boundary_normal_form.tsv'), ...
    'FileType','text','Delimiter','\t');

rng(cfg.RandomSeed);
directions = {r,-r,primary.TransientAudit.Physical.OptimalInput, ...
    -primary.TransientAudit.Physical.OptimalInput};
directionNames = {'critical+','critical-','transient+','transient-'};
for randomIndex = 1:3
    random = randn(size(fixed)); random=random/norm(random);
    directions{end+1}=random; %#ok<AGROW>
    directionNames{end+1}=sprintf('random%d',randomIndex); %#ok<AGROW>
end
gainTest = [gammaCritical-0.02,gammaCritical+0.02];
time = (0:0.5:200)';
trajectoryRows = {};
for gainIndex = 1:2
    testPhi = @(x)l6ns_phi(x,0,context,[1 1],gainTest(gainIndex));
    for directionIndex = 1:numel(directions)
        direction = directions{directionIndex}/norm(directions{directionIndex});
        state = fixed+1e-4*stateScale*direction;
        initialDeviation = norm(state-fixed);
        maximumDeviation = initialDeviation;
        for ti = 2:numel(time)
            dt=time(ti)-time(ti-1);
            state=state+(dt/cfg.TauMs)*(testPhi(state)-state);
            maximumDeviation=max(maximumDeviation,norm(state-fixed));
        end
        residual=norm(testPhi(state)-state)/max(norm(state),eps);
        bounded=all(isfinite(state)) && max(abs(state))<1e4;
        n=data.PopulationSize;
        trajectoryRows(end+1,:)={gainTest(gainIndex),directionNames{directionIndex}, ...
            norm(state-fixed)/initialDeviation,maximumDeviation,residual,bounded, ...
            mean(state(1:n)),mean(state(n+(1:n))),mean(state(2*n+(1:n)))}; %#ok<AGROW>
    end
end
trajectories = cell2table(trajectoryRows,'VariableNames', ...
    {'gammaI','direction','finalGrowth','maximumDeviation','finalResidual', ...
    'bounded','finalMeanS','finalMeanC','finalMeanI'});
writetable(trajectories,fullfile(outputDir,'inhibition_boundary_directions.tsv'), ...
    'FileType','text','Delimiter','\t');

summary = table(gammaCritical,real(lambda),imag(lambda),separation,real(beta), ...
    median(quadratic(end-2:end)),median(cubic(end-2:end)), ...
    abs(imag(lambda))<1e-7 && separation>1e-5, ...
    norm(phi(fixed)-fixed)/max(norm(fixed),eps), ...
    'VariableNames',{'gammaICritical','lambdaReal','lambdaImag', ...
    'spectralSeparation','parameterSensitivity','quadraticCoefficient', ...
    'cubicCoefficient','simpleRealCrossing','relativeFixedPointResidual'});
writetable(summary,fullfile(outputDir,'inhibition_boundary_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
boundary = struct('Status','boundary_found','Summary',summary, ...
    'NormalForm',normalForm,'Trajectories',trajectories);
end
