function result = l6ns_two_pathway_plan(cfg,data,context)
% Execute the proposal's canonical fixed-point-preserving two-pathway test.

outputDir = fullfile(cfg.OutputRoot,'two_pathway_plan');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
pathway = l6ns_two_pathway_data(data);

fixed = context.FixedPoint(:);
rawContext = context;
rawContext.FixedPointCorrection = [];
rawBaseline = l6ns_phi(fixed,0,rawContext,[1 1],1);
context.FixedPointCorrection = fixed-rawBaseline;
operators = l6ns_control_observation(context,cfg);

audit = local_audit(cfg,pathway,context,fixed,outputDir);
weights = local_population_weights(fixed,pathway.PopulationSize);

l6Branch = l6ns_pathway_branch(pathway,cfg.Gamma6Grid,'L6',cfg);
iBranch = l6ns_pathway_branch(pathway,cfg.GammaIGrid,'I',cfg);
writetable(l6Branch.Table,fullfile(outputDir,'l6_branch.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(iBranch.Table,fullfile(outputDir,'inhibition_branch.tsv'), ...
    'FileType','text','Delimiter','\t');

transientAudit = local_transient_audit( ...
    pathway.JBaseline,cfg,weights,outputDir);
matched = local_matched_stability( ...
    cfg,pathway,l6Branch,iBranch,operators,weights,outputDir);
interaction = local_interaction_grid( ...
    cfg,pathway,l6Branch,iBranch,operators,weights,outputDir);
loop = local_ordered_loop(cfg,pathway,l6Branch,outputDir);
synthetic = local_synthetic_readout(matched,operators,outputDir);
opposite = local_opposite_effect_summary( ...
    cfg,pathway,l6Branch,iBranch,matched,interaction,outputDir);
inhibitionFate = local_inhibition_fate( ...
    cfg,pathway,iBranch,context,fixed,outputDir);

local_plot_branches(cfg,l6Branch,iBranch,outputDir);
local_plot_matched(matched,outputDir);
local_plot_grid(interaction,outputDir);
local_plot_mode_maps(pathway,l6Branch,iBranch,data,outputDir);
local_plot_matched_singular_maps(matched,data,outputDir);

result = struct('Config',cfg,'Pathway',pathway,'Audit',audit, ...
    'Operators',operators,'L6Branch',l6Branch,'InhibitionBranch',iBranch, ...
    'TransientAudit',transientAudit,'Matched',matched, ...
    'Interaction',interaction,'OrderedLoop',loop,'Synthetic',synthetic);
result.OppositeEffect = opposite;
result.InhibitionFate = inhibitionFate;
save(fullfile(outputDir,'two_pathway_result.mat'),'result','-v7.3');
end

function audit = local_audit(cfg,pathway,context,fixed,outputDir)
gainPairs = [1 1;0 1;1 0;0.8 1.2];
rng(cfg.RandomSeed);
directionCount = 3;
jvpRows = {};
fixedRows = cell(size(gainPairs,1),1);
step = 2e-6*max(1,norm(fixed)/sqrt(numel(fixed)));
for pairIndex = 1:size(gainPairs,1)
    gamma6 = gainPairs(pairIndex,1);
    gammaI = gainPairs(pairIndex,2);
    jacobian = l6ns_pathway_jacobian(pathway,gamma6,gammaI);
    phi = @(x)l6ns_phi(x,1-gamma6,context,[1 1],gammaI);
    fixedResidual = norm(phi(fixed)-fixed)/max(norm(fixed),eps);
    fixedRows{pairIndex} = {gamma6,gammaI,fixedResidual};
    for directionIndex = 1:directionCount
        z = randn(size(fixed)); z = z/norm(z);
        numeric = (phi(fixed+step*z)-phi(fixed-step*z))/(2*step);
        analytic = jacobian*z;
        jvpError = norm(numeric-analytic)/max(norm(analytic),eps);
        affineAction = pathway.JRest*z+gamma6*(pathway.J6*z)+gammaI*(pathway.JI*z);
        decompositionError = norm(analytic-affineAction)/max(norm(analytic),eps);
        jvpRows(end+1,:) = {gamma6,gammaI,directionIndex,jvpError, ...
            decompositionError}; %#ok<AGROW>
    end
end
jvpTable = cell2table(jvpRows,'VariableNames', ...
    {'gamma6','gammaI','direction','relativeJvpError','decompositionJvpError'});
fixedTable = cell2table(vertcat(fixedRows{:}),'VariableNames', ...
    {'gamma6','gammaI','relativeFixedPointResidual'});
writetable(jvpTable,fullfile(outputDir,'pathway_jvp_audit.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(fixedTable,fullfile(outputDir,'fixed_point_preservation.tsv'), ...
    'FileType','text','Delimiter','\t');

n = pathway.PopulationSize;
blockRows = {};
labels = {'S','C','I'};
for target = 1:3
    rows = (target-1)*n+(1:n);
    block = pathway.JI(rows,pathway.IColumns);
    values = nonzeros(block);
    blockRows(end+1,:) = {labels{target},norm(block,'fro'), ...
        mean(values>0),mean(values<0),min(values),max(values)}; %#ok<AGROW>
end
blockTable = cell2table(blockRows,'VariableNames', ...
    {'target','frobeniusNorm','positiveEntryFraction','negativeEntryFraction', ...
    'minimumEntry','maximumEntry'});
writetable(blockTable,fullfile(outputDir,'inhibitory_block_sign_audit.tsv'), ...
    'FileType','text','Delimiter','\t');

summary = table(cfg.TauMs,pathway.ReconstructionError, ...
    pathway.J6ISourceFraction,pathway.JINonISourceFraction, ...
    max(jvpTable.relativeJvpError),max(fixedTable.relativeFixedPointResidual), ...
    'VariableNames',{'tauMs','pathwayReconstructionError','J6ISourceFraction', ...
    'JINonISourceFraction','maximumJvpError','maximumFixedPointResidual'});
writetable(summary,fullfile(outputDir,'operator_audit_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
audit = struct('Summary',summary,'JVP',jvpTable,'FixedPoint',fixedTable, ...
    'InhibitoryBlocks',blockTable);
end

function transient = local_transient_audit(jacobian,cfg,weights,outputDir)
stepGrid = (0:20)';
oldGain = nan(size(stepGrid));
for k = 1:numel(stepGrid)
    oldGain(k) = local_discrete_gain(jacobian,stepGrid(k),weights,cfg.RandomSeed);
end
physical = l6ns_continuous_transient( ...
    jacobian,cfg.PhysicalTimeGridMs,cfg,weights,cfg.TauMs);
physicalEuclidean = l6ns_continuous_transient( ...
    jacobian,cfg.PhysicalTimeGridMs,cfg,[],cfg.TauMs);
oldTable = table(stepGrid,oldGain,'VariableNames',{'iterationStep','sigma1JPower'});
physicalTable = table(physical.TimeGrid,physical.Gain, ...
    'VariableNames',{'timeMs','sigma1ExpLt'});
writetable(oldTable,fullfile(outputDir,'old_iteration_transient.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(physicalTable,fullfile(outputDir,'physical_transient.tsv'), ...
    'FileType','text','Delimiter','\t');
summary = table(max(oldGain),stepGrid(find(oldGain==max(oldGain),1)), ...
    physical.PeakGain,physical.PeakTime,physical.ForwardSingularResidual, ...
    physical.AdjointSingularResidual,physicalEuclidean.PeakGain, ...
    physicalEuclidean.PeakTime,physicalEuclidean.ForwardSingularResidual, ...
    physicalEuclidean.AdjointSingularResidual, ...
    physical.PeakGain/physicalEuclidean.PeakGain,'VariableNames', ...
    {'oldPeakGain','oldPeakStep','physicalPeakGain','physicalPeakTimeMs', ...
    'forwardSingularResidual','adjointSingularResidual','euclideanPeakGain', ...
    'euclideanPeakTimeMs','euclideanForwardResidual','euclideanAdjointResidual', ...
    'scaledToEuclideanPeakGainRatio'});
writetable(summary,fullfile(outputDir,'transient_audit_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[100 100 1100 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(stepGrid,oldGain,'-o','LineWidth',1.4); grid on;
xlabel('iteration step n'); ylabel('\sigma_1(J^n)'); title('Old discrete metric');
nexttile; plot(physical.TimeGrid,physical.Gain,'-o','LineWidth',1.4); hold on;
plot(physicalEuclidean.TimeGrid,physicalEuclidean.Gain,'--','LineWidth',1.3); grid on;
xlabel('physical time (ms)'); ylabel('\sigma_1(e^{Lt})');
title(sprintf('Physical metric, \\tau=%.4f ms',cfg.TauMs));
legend({'population RMS scaled','Euclidean'},'Location','best');
l6ns_save_figure(fig,outputDir,'figure3_old_vs_physical_transient'); close(fig);
transient = struct('Old',oldTable,'Physical',physical, ...
    'PhysicalEuclidean',physicalEuclidean,'Summary',summary);
end

function gain = local_discrete_gain(jacobian,power,weights,seedValue)
n = size(jacobian,1);
w = spdiags(weights,0,n,n);
iw = spdiags(1./weights,0,n,n);
matrix = w*jacobian*iw;
rng(seedValue);
v = randn(n,1); v = v/norm(v);
for iteration = 1:10
    u = v;
    for k = 1:power; u = matrix*u; end
    u = u/max(norm(u),eps);
    vNew = u;
    for k = 1:power; vNew = matrix'*vNew; end
    vNew = vNew/max(norm(vNew),eps);
    if abs(vNew'*v)>1-1e-8; v=vNew; break; end
    v = vNew;
end
u = v;
for k = 1:power; u = matrix*u; end
gain = norm(u);
end

function matched = local_matched_stability(cfg,pathway,l6Branch,iBranch,operators,weights,outputDir)
l6Stable = l6Branch.Table.spectralAbscissaJ<1;
iStable = iBranch.Table.spectralAbscissaJ<1;
low = max(min(l6Branch.Table.spectralAbscissaJ(l6Stable)), ...
    min(iBranch.Table.spectralAbscissaJ(iStable)));
high = min(max(l6Branch.Table.spectralAbscissaJ(l6Stable)), ...
    max(iBranch.Table.spectralAbscissaJ(iStable)));
if ~(isfinite(low)&&isfinite(high)&&high>low)
    error('L6 and inhibition sweeps have no overlapping stable alpha range.');
end
fraction = linspace(0.2,0.85,cfg.MatchedStabilityCount);
targets = low+(high-low)*fraction;
rows = cell(numel(targets),1);
pairs = cell(numel(targets),1);
for targetIndex = 1:numel(targets)
    target = targets(targetIndex);
    gamma6 = local_match_gain(pathway,l6Branch.Table.gamma, ...
        l6Branch.Table.spectralAbscissaJ,target,'L6',cfg);
    gammaI = local_match_gain(pathway,iBranch.Table.gamma, ...
        iBranch.Table.spectralAbscissaJ,target,'I',cfg);
    j6 = l6ns_pathway_jacobian(pathway,gamma6,1);
    jI = l6ns_pathway_jacobian(pathway,1,gammaI);
    modes6 = l6ns_eigenpairs(j6,6,'largestreal',cfg);
    modesI = l6ns_eigenpairs(jI,6,'largestreal',cfg);
    alpha6 = max(real(modes6.Lambda));
    alphaI = max(real(modesI.Lambda));
    physical6 = l6ns_continuous_transient(j6,cfg.PhysicalTimeGridMs, ...
        cfg,weights,cfg.TauMs);
    physicalI = l6ns_continuous_transient(jI,cfg.PhysicalTimeGridMs, ...
        cfg,weights,cfg.TauMs);
    constrained6 = l6ns_constrained_transient(j6,cfg.PhysicalTimeGridMs,cfg,operators);
    constrainedI = l6ns_constrained_transient(jI,cfg.PhysicalTimeGridMs,cfg,operators);
    qR6 = orth(modes6.Right(:,1:4)); qRI = orth(modesI.Right(:,1:4));
    qL6 = orth(modes6.Left(:,1:4)); qLI = orth(modesI.Left(:,1:4));
    rightCosine = min(svd(qR6'*qRI));
    leftCosine = min(svd(qL6'*qLI));
    access6 = l6ns_mode_accessibility(modes6,operators,cfg);
    accessI = l6ns_mode_accessibility(modesI,operators,cfg);
    kernelDifference = norm(constrained6.Kernels(:)-constrainedI.Kernels(:)) / ...
        max(norm(constrained6.Kernels(:)),eps);
    rows{targetIndex} = {target,gamma6,gammaI,alpha6,alphaI,abs(alpha6-alphaI), ...
        cfg.TauMs/(1-alpha6),cfg.TauMs/(1-alphaI), ...
        modes6.ConditionNumber(1),modesI.ConditionNumber(1),rightCosine,leftCosine, ...
        physical6.PeakGain,physicalI.PeakGain,physical6.PeakTime,physicalI.PeakTime, ...
        constrained6.PeakGain,constrainedI.PeakGain, ...
        constrained6.PeakTimeMs,constrainedI.PeakTimeMs, ...
        constrained6.PeakAfterZeroGain,constrainedI.PeakAfterZeroGain, ...
        constrained6.PeakAfterZeroTimeMs,constrainedI.PeakAfterZeroTimeMs, ...
        kernelDifference, ...
        access6.labIdentifiability(1),accessI.labIdentifiability(1)};
    pairs{targetIndex} = struct('Target',target,'Gamma6',gamma6,'GammaI',gammaI, ...
        'Modes6',modes6,'ModesI',modesI, ...
        'Physical6',physical6,'PhysicalI',physicalI, ...
        'Constrained6',constrained6,'ConstrainedI',constrainedI, ...
        'Accessibility6',access6,'AccessibilityI',accessI);
end
summary = cell2table(vertcat(rows{:}),'VariableNames', ...
    {'targetLambdaReal','gamma6','gammaI','lambda6Real','lambdaIReal', ...
    'matchingError','recovery6Ms','recoveryIMs','condition6','conditionI', ...
    'minimumRightSubspaceCosine','minimumLeftSubspaceCosine', ...
    'physicalGain6','physicalGainI','physicalPeak6Ms','physicalPeakIMs', ...
    'constrainedGain6','constrainedGainI','constrainedPeak6Ms', ...
    'constrainedPeakIMs','constrainedAfterZeroGain6','constrainedAfterZeroGainI', ...
    'constrainedAfterZeroPeak6Ms','constrainedAfterZeroPeakIMs', ...
    'relativeKernelDifference','labScore6','labScoreI'});
writetable(summary,fullfile(outputDir,'matched_stability_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
matched = struct('Summary',summary,'Pairs',{pairs});
end

function gamma = local_match_gain(pathway,gammaGrid,alphaGrid,target,axisName,cfg)
crossing = find((alphaGrid(1:end-1)-target).*(alphaGrid(2:end)-target)<=0);
if isempty(crossing)
    [~,nearest] = min(abs(alphaGrid-target));
    gamma = gammaGrid(nearest);
    return
end
[~,choice] = min(abs((gammaGrid(crossing)+gammaGrid(crossing+1))/2-1));
index = crossing(choice);
low = gammaGrid(index); high = gammaGrid(index+1);
fLow = alphaGrid(index)-target;
for iteration = 1:24
    middle = (low+high)/2;
    if strcmp(axisName,'L6')
        matrix = l6ns_pathway_jacobian(pathway,middle,1);
    else
        matrix = l6ns_pathway_jacobian(pathway,1,middle);
    end
    fMiddle = l6ns_max_real(matrix,cfg)-target;
    if fLow*fMiddle<=0
        high=middle;
    else
        low=middle; fLow=fMiddle;
    end
end
gamma = (low+high)/2;
end

function interaction = local_interaction_grid(cfg,pathway,l6Branch,iBranch,operators,weights,outputDir)
[~,l6BaseIndex] = min(abs(l6Branch.Table.gamma-1));
[~,iBaseIndex] = min(abs(iBranch.Table.gamma-1));
s6 = l6Branch.Table.sensitivity6Real(l6BaseIndex);
sI = iBranch.Table.sensitivityIReal(iBaseIndex);
slope = -s6/sI;
gamma6 = 1+cfg.InteractionGamma6Offsets;
iHalfWidth = min(0.3,max(0.08,abs(slope)*max(abs(cfg.InteractionGamma6Offsets))));
gammaI = 1+linspace(-iHalfWidth,iHalfWidth,5);

shape = [numel(gammaI),numel(gamma6)];
alpha = nan(shape); numerical = nan(shape); condition = nan(shape);
physicalGain = nan(shape); constrainedGain = nan(shape); peakTime = nan(shape);
for iIndex = 1:numel(gammaI)
    for sixIndex = 1:numel(gamma6)
        matrix = l6ns_pathway_jacobian(pathway,gamma6(sixIndex),gammaI(iIndex));
        modes = l6ns_eigenpairs(matrix,4,'largestreal',cfg);
        alpha(iIndex,sixIndex) = max(real(modes.Lambda));
        condition(iIndex,sixIndex) = modes.ConditionNumber(1);
        generator = (matrix-speye(size(matrix,1)))/cfg.TauMs;
        numerical(iIndex,sixIndex) = real(eigs((generator+generator')/2,1,'largestreal'));
        gridCfg = cfg;
        gridCfg.TransientPowerIterations = min(10,cfg.TransientPowerIterations);
        gridCfg.TransientRefinePeak = false;
        physical = l6ns_continuous_transient(matrix,cfg.GridTimeGridMs, ...
            gridCfg,weights,cfg.TauMs);
        constrained = l6ns_constrained_transient(matrix,cfg.GridTimeGridMs,cfg,operators);
        physicalGain(iIndex,sixIndex) = physical.PeakGain;
        constrainedGain(iIndex,sixIndex) = constrained.PeakAfterZeroGain;
        peakTime(iIndex,sixIndex) = physical.PeakTime;
    end
end

baselineRow = find(abs(gammaI-1)<1e-12,1);
baselineColumn = find(abs(gamma6-1)<1e-12,1);
if isempty(baselineRow) || isempty(baselineColumn)
    error('Interaction grid must include (1,1).');
end
deltaAlpha = alpha-alpha(baselineRow,:)-alpha(:,baselineColumn)+alpha(baselineRow,baselineColumn);
deltaPhysical = physicalGain-physicalGain(baselineRow,:)- ...
    physicalGain(:,baselineColumn)+physicalGain(baselineRow,baselineColumn);
deltaConstrained = constrainedGain-constrainedGain(baselineRow,:)- ...
    constrainedGain(:,baselineColumn)+constrainedGain(baselineRow,baselineColumn);
deltaCondition = condition-condition(baselineRow,:)- ...
    condition(:,baselineColumn)+condition(baselineRow,baselineColumn);

[g6Mesh,gIMesh] = meshgrid(gamma6,gammaI);
gridTable = table(g6Mesh(:),gIMesh(:),alpha(:),(1-alpha(:)),numerical(:), ...
    condition(:),physicalGain(:),constrainedGain(:),peakTime(:), ...
    deltaAlpha(:),deltaPhysical(:),deltaConstrained(:),deltaCondition(:), ...
    'VariableNames',{'gamma6','gammaI','lambdaMaxReal','stabilityMarginJ', ...
    'numericalAbscissaPerMs','criticalConditionNumber','physicalPeakGain', ...
    'constrainedPeakGain','physicalPeakTimeMs','nonadditiveAlpha', ...
    'nonadditivePhysicalGain','nonadditiveConstrainedGain','nonadditiveCondition'});
writetable(gridTable,fullfile(outputDir,'interaction_grid.tsv'), ...
    'FileType','text','Delimiter','\t');
compensation = table(s6,sI,slope,iHalfWidth,'VariableNames', ...
    {'s6','sI','predictedDGammaIDGamma6','gammaIHalfWidth'});
writetable(compensation,fullfile(outputDir,'compensation_direction.tsv'), ...
    'FileType','text','Delimiter','\t');
interaction = struct('Gamma6',gamma6,'GammaI',gammaI,'Alpha',alpha, ...
    'NumericalAbscissa',numerical,'Condition',condition, ...
    'PhysicalGain',physicalGain,'ConstrainedGain',constrainedGain, ...
    'PeakTimeMs',peakTime,'DeltaAlpha',deltaAlpha, ...
    'DeltaPhysicalGain',deltaPhysical,'DeltaConstrainedGain',deltaConstrained, ...
    'DeltaCondition',deltaCondition,'Table',gridTable,'Compensation',compensation);
end

function loop = local_ordered_loop(cfg,pathway,l6Branch,outputDir)
[~,baseIndex] = min(abs(l6Branch.Table.gamma-1));
baseline = l6ns_pathway_jacobian(pathway,1,1);
modes = l6ns_eigenpairs(baseline,8,'largestreal',cfg);
q = orth(modes.Right);
iThen6 = q'*(pathway.JI*(pathway.J6*q));
sixThenI = q'*(pathway.J6*(pathway.JI*q));
commutator = iThen6-sixThenI;
fullScale = max(norm(iThen6,'fro')+norm(sixThenI,'fro'),eps);
summary = table(size(q,2),norm(iThen6,'fro'),norm(sixThenI,'fro'), ...
    norm(commutator,'fro'),norm(commutator,'fro')/fullScale, ...
    l6Branch.Table.sensitivity6Real(baseIndex), ...
    'VariableNames',{'subspaceDimension','normJIJ6','normJ6JI', ...
    'commutatorNorm','relativeOrderedLoopMismatch','baselineS6'});
writetable(summary,fullfile(outputDir,'ordered_pathway_loop.tsv'), ...
    'FileType','text','Delimiter','\t');
loop = struct('Summary',summary,'Basis',q,'JIJ6',iThen6, ...
    'J6JI',sixThenI,'Commutator',commutator);
end

function synthetic = local_synthetic_readout(matched,operators,outputDir)
pair = matched.Pairs{end};
time = pair.Constrained6.TimeGridMs;
k6 = pair.Constrained6.Kernels;
kI = pair.ConstrainedI.Kernels;
rows = {};
for observation = 1:numel(operators.ObservationNames)
    for control = 1:numel(operators.ControlNames)
        for ti = 1:numel(time)
            rows(end+1,:) = {operators.ObservationNames{observation}, ...
                operators.ControlNames{control},time(ti), ...
                k6(observation,control,ti),kI(observation,control,ti)}; %#ok<AGROW>
        end
    end
end
kernelTable = cell2table(rows,'VariableNames', ...
    {'observation','control','timeMs','l6Kernel','inhibitionKernel'});
writetable(kernelTable,fullfile(outputDir,'synthetic_response_kernels.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[80 80 1350 820]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
controlIndex = [1 3 5];
for row = 1:2
    observation = row;
    for column = 1:3
        nexttile;
        plot(time,squeeze(k6(observation,controlIndex(column),:)),'-','LineWidth',1.6); hold on;
        plot(time,squeeze(kI(observation,controlIndex(column),:)),'--','LineWidth',1.6);
        grid on; xlabel('time (ms)'); ylabel('response');
        title(sprintf('%s, %s',operators.ObservationNames{observation}, ...
            operators.ControlNames{controlIndex(column)}),'Interpreter','none');
        if row==1 && column==1; legend({'L6 matched','I matched'},'Location','best'); end
    end
end
l6ns_save_figure(fig,outputDir,'figure4_matched_response_kernels'); close(fig);
synthetic = struct('KernelTable',kernelTable,'PairIndex',numel(matched.Pairs));
end

function summary = local_opposite_effect_summary(cfg,pathway,l6Branch,iBranch,matched,interaction,outputDir)
[~,l6Index] = min(abs(l6Branch.Table.gamma-1));
[~,iIndex] = min(abs(iBranch.Table.gamma-1));
baseline = l6ns_max_real(l6ns_pathway_jacobian(pathway,1,1),cfg);
frozenL6 = l6ns_max_real(l6ns_pathway_jacobian(pathway,0,1),cfg);
frozenI = l6ns_max_real(l6ns_pathway_jacobian(pathway,1,0),cfg);
s6 = l6Branch.Table.sensitivity6Real(l6Index);
sI = iBranch.Table.sensitivityIReal(iIndex);
summary = table(s6,sI,sign(s6)*sign(sI)<0,baseline,frozenL6,frozenI, ...
    frozenL6-baseline,frozenI-baseline, ...
    interaction.Compensation.predictedDGammaIDGamma6, ...
    max(matched.Summary.matchingError), ...
    mean(matched.Summary.physicalGain6-matched.Summary.physicalGainI), ...
    mean(matched.Summary.constrainedAfterZeroGain6- ...
        matched.Summary.constrainedAfterZeroGainI), ...
    'VariableNames',{'s6','sI','oppositeLocalSpectralSigns','baselineLambdaReal', ...
    'frozenL6LambdaReal','frozenILambdaReal','freezeL6DeltaLambda', ...
    'freezeIDeltaLambda','constantStabilitySlope','maximumMatchingError', ...
    'meanMatchedPhysicalGainDifference','meanMatchedConstrainedGainDifference'});
writetable(summary,fullfile(outputDir,'opposite_effect_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
end

function fate = local_inhibition_fate(cfg,pathway,iBranch,context,fixed,outputDir)
alpha = iBranch.Table.spectralAbscissaJ;
gamma = iBranch.Table.gamma;
bracket = find((alpha(1:end-1)-1).*(alpha(2:end)-1)<=0,1,'last');
if isempty(bracket)
    fate = struct('Status','no_boundary_in_sweep');
    return
end
low = gamma(bracket); high = gamma(bracket+1);
fLow = alpha(bracket)-1;
for iteration = 1:28
    middle = (low+high)/2;
    fMiddle = l6ns_max_real(l6ns_pathway_jacobian(pathway,1,middle),cfg)-1;
    if fLow*fMiddle<=0
        high=middle;
    else
        low=middle; fLow=fMiddle;
    end
end
gammaCritical = (low+high)/2;
matrixCritical = l6ns_pathway_jacobian(pathway,1,gammaCritical);
modes = l6ns_eigenpairs(matrixCritical,6,'largestreal',cfg);
criticalLambda = modes.Lambda(1);
r = real(modes.Right(:,1)); r=r/norm(r);
epsilon = 1e-4*max(1,norm(fixed));
gammaTest = [gammaCritical-0.02,gammaCritical+0.02];
signTest = [-1 1];
timeMs = (0:0.5:100)';
rows = {};
trajectories = cell(2,2);
for gammaIndex = 1:2
    gammaI = gammaTest(gammaIndex);
    phi = @(x)l6ns_phi(x,0,context,[1 1],gammaI);
    for signIndex = 1:2
        state = fixed+signTest(signIndex)*epsilon*r;
        deviation = nan(numel(timeMs),1);
        deviation(1)=norm(state-fixed);
        for ti = 2:numel(timeMs)
            dt = timeMs(ti)-timeMs(ti-1);
            state = state+(dt/cfg.TauMs)*(phi(state)-state);
            deviation(ti)=norm(state-fixed);
        end
        finalResidual = norm(phi(state)-state)/max(norm(state),eps);
        growth = deviation(end)/max(deviation(1),eps);
        bounded = all(isfinite(state)) && max(abs(state))<1e4;
        rows(end+1,:) = {gammaI,signTest(signIndex),growth,max(deviation), ...
            finalResidual,bounded}; %#ok<AGROW>
        trajectories{gammaIndex,signIndex}=deviation;
    end
end
summary = cell2table(rows,'VariableNames', ...
    {'gammaI','perturbationSign','finalGrowth','maximumDeviation', ...
    'finalFixedPointResidual','bounded'});
writetable(summary,fullfile(outputDir,'inhibition_first_boundary_fate.tsv'), ...
    'FileType','text','Delimiter','\t');
boundary = table(gammaCritical,real(criticalLambda),imag(criticalLambda), ...
    modes.ConditionNumber(1),'VariableNames', ...
    {'gammaICritical','lambdaReal','lambdaImag','conditionNumber'});
writetable(boundary,fullfile(outputDir,'inhibition_first_boundary.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[100 100 900 500]);
legendLabels = {'below boundary (-)','below boundary (+)', ...
    'above boundary (-)','above boundary (+)'};
lineIndex = 0;
for gammaIndex = 1:2
    for signIndex = 1:2
        lineIndex = lineIndex+1;
        semilogy(timeMs,trajectories{gammaIndex,signIndex},'LineWidth',1.4, ...
            'DisplayName',legendLabels{lineIndex}); hold on;
    end
end
grid on; xlabel('physical time (ms)'); ylabel('||f(t)-f_*||_2');
legend('show','Location','best');
title(sprintf('First inhibition boundary \\gamma_I^c=%.6f',gammaCritical));
l6ns_save_figure(fig,outputDir,'figure6_inhibition_first_boundary_fate'); close(fig);
fate = struct('Status','boundary_found','Boundary',boundary,'Summary',summary, ...
    'TimeMs',timeMs,'Trajectories',{trajectories});
end

function weights = local_population_weights(fixed,n)
scale = [norm(fixed(1:n))/sqrt(n),norm(fixed(n+(1:n)))/sqrt(n), ...
    norm(fixed(2*n+(1:n)))/sqrt(n)];
scale = max(scale,1e-8);
weights = [ones(n,1)/scale(1);ones(n,1)/scale(2);ones(n,1)/scale(3)];
end

function local_plot_branches(cfg,l6Branch,iBranch,outputDir)
fig = figure('Visible','off','Color','w','Position',[100 100 1350 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(l6Branch.Table.gamma,l6Branch.Table.spectralAbscissaJ,'-o','LineWidth',1.4); hold on;
plot(l6Branch.Table.gamma,l6Branch.Table.lambdaReal,'--','LineWidth',1.1);
yline(1,'r--'); grid on; xlabel('\gamma_6'); ylabel('max Re \lambda(J)'); title('L6 branch');
nexttile; plot(iBranch.Table.gamma,iBranch.Table.spectralAbscissaJ,'-o','LineWidth',1.4); hold on;
plot(iBranch.Table.gamma,iBranch.Table.lambdaReal,'--','LineWidth',1.1);
yline(1,'r--'); grid on; xlabel('\gamma_I'); ylabel('max Re \lambda(J)'); title('Inhibition branch');
nexttile; plot(l6Branch.Table.gamma,l6Branch.Table.L6toSReal,'-o'); hold on;
plot(l6Branch.Table.gamma,l6Branch.Table.L6toCReal,'-o');
plot(l6Branch.Table.gamma,l6Branch.Table.L6toIReal,'-o'); grid on;
xlabel('\gamma_6'); ylabel('Re(y^*J_{6\rightarrow p}x)'); legend({'S','C','I'});
nexttile; plot(iBranch.Table.gamma,iBranch.Table.ItoSReal,'-o'); hold on;
plot(iBranch.Table.gamma,iBranch.Table.ItoCReal,'-o');
plot(iBranch.Table.gamma,iBranch.Table.ItoIReal,'-o'); grid on;
xlabel('\gamma_I'); ylabel('Re(y^*J_{I\rightarrow p}x)'); legend({'S','C','I'});
sgtitle(sprintf('Physical generator L=(J-I)/%.4f ms',cfg.TauMs));
l6ns_save_figure(fig,outputDir,'figure2_l6_and_inhibition_branches'); close(fig);
end

function local_plot_matched(matched,outputDir)
t = matched.Summary;
fig = figure('Visible','off','Color','w','Position',[100 100 1350 430]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; plot(t.targetLambdaReal,t.physicalGain6,'-o','LineWidth',1.4); hold on;
plot(t.targetLambdaReal,t.physicalGainI,'--o','LineWidth',1.4); grid on;
xlabel('matched max Re \lambda'); ylabel('physical peak gain'); legend({'L6','I'});
nexttile; plot(t.targetLambdaReal,t.constrainedAfterZeroGain6,'-o','LineWidth',1.4); hold on;
plot(t.targetLambdaReal,t.constrainedAfterZeroGainI,'--o','LineWidth',1.4); grid on;
xlabel('matched max Re \lambda'); ylabel('constrained gain, t>0');
nexttile; semilogy(t.targetLambdaReal,t.condition6,'-o','LineWidth',1.4); hold on;
semilogy(t.targetLambdaReal,t.conditionI,'--o','LineWidth',1.4); grid on;
xlabel('matched max Re \lambda'); ylabel('critical condition number');
l6ns_save_figure(fig,outputDir,'figure4_matched_stability_phenotypes'); close(fig);
end

function local_plot_grid(interaction,outputDir)
fig = figure('Visible','off','Color','w','Position',[80 80 1400 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
values = {interaction.Alpha,interaction.PhysicalGain, ...
    interaction.ConstrainedGain,interaction.Condition};
titles = {'max Re \lambda(J)','physical peak gain','constrained peak gain', ...
    'critical condition number'};
for index = 1:4
    nexttile; imagesc(interaction.Gamma6,interaction.GammaI,values{index});
    axis xy; colorbar; xlabel('\gamma_6'); ylabel('\gamma_I'); title(titles{index});
end
l6ns_save_figure(fig,outputDir,'figure4_interaction_grid'); close(fig);
end

function local_plot_mode_maps(pathway,l6Branch,iBranch,data,outputDir)
[~,i6] = min(abs(l6Branch.Table.gamma-1));
[~,iI] = min(abs(iBranch.Table.gamma-1));
vectors = {l6Branch.Right{i6},l6Branch.Left{i6},iBranch.Right{iI},iBranch.Left{iI}};
rowLabels = {'L6 right','L6 left','I right','I left'};
n = pathway.PopulationSize;
fig = figure('Visible','off','Color','w','Position',[50 50 1450 1050]);
tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
populationLabels = {'S','C','I','E'};
for row = 1:4
    v = real(vectors{row});
    maps = {reshape(v(1:n),data.MapSize),reshape(v(n+(1:n)),data.MapSize), ...
        reshape(v(2*n+(1:n)),data.MapSize)};
    maps{4} = (1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit = max(cellfun(@(x)max(abs(x(:))),maps));
    for column = 1:4
        nexttile; imagesc(maps{column}); axis image off; colorbar;
        if limit>0; clim([-limit limit]); end
        if row==1; title(populationLabels{column}); end
        if column==1
            text(-0.12,0.5,rowLabels{row},'Units','normalized', ...
                'Rotation',90,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold');
        end
    end
end
l6ns_save_figure(fig,outputDir,'figure4_l6_inhibition_left_right_maps'); close(fig);
end

function local_plot_matched_singular_maps(matched,data,outputDir)
pair = matched.Pairs{end};
vectors = {pair.Physical6.OptimalInput,pair.Physical6.PeakOutput, ...
    pair.PhysicalI.OptimalInput,pair.PhysicalI.PeakOutput};
rowLabels = {'L6 input','L6 output','I input','I output'};
populationLabels = {'S','C','I','E'};
n = data.PopulationSize;
fig = figure('Visible','off','Color','w','Position',[50 50 1450 1050]);
tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
for row = 1:4
    v = real(vectors{row});
    maps = {reshape(v(1:n),data.MapSize),reshape(v(n+(1:n)),data.MapSize), ...
        reshape(v(2*n+(1:n)),data.MapSize)};
    maps{4} = (1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit = max(cellfun(@(x)max(abs(x(:))),maps));
    for column = 1:4
        nexttile; imagesc(maps{column}); axis image off; colorbar;
        if limit>0; clim([-limit limit]); end
        if row==1; title(populationLabels{column}); end
        if column==1
            text(-0.12,0.5,rowLabels{row},'Units','normalized', ...
                'Rotation',90,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontWeight','bold');
        end
    end
end
l6ns_save_figure(fig,outputDir,'figure3_matched_singular_input_output_maps'); close(fig);
end
