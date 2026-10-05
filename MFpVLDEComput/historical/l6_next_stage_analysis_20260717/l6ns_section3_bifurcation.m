function result = l6ns_section3_bifurcation(cfg, data, section1, context)
% Section 3: classify the lambda=1 event of dot(f)=Phi(f)-f.

sectionDir = fullfile(cfg.OutputRoot, 'section3_bifurcation');
if ~exist(sectionDir, 'dir'); mkdir(sectionDir); end

wc = section1.Summary.directWc;
alphaC = 1 - wc;
jCritical = data.A + alphaC * data.B;
modes = l6ns_eigenpairs(jCritical, 16, 'largestreal', cfg);
[~,criticalIndex] = min(abs(modes.Lambda - 1));
lambdaCritical = modes.Lambda(criticalIndex);
r = modes.Right(:,criticalIndex);
l = modes.Left(:,criticalIndex);
r = real(r / norm(r));
l = real(l / conj(l'*r));
l = l / (l'*r);

other = modes.Lambda;
other(criticalIndex) = [];
spectralSeparation = min(abs(other - lambdaCritical));
secondLeadingReal = max(real(other));
beta = -l' * data.B * r;

fixed = context.FixedPoint;
fixedPointDifference = norm(fixed-data.FixedPoint) / max(norm(fixed),eps);
rawContext = context;
rawFixedResponse = l6ns_phi(fixed,1,rawContext);
rawFixedResidual = norm(rawFixedResponse-fixed) / max(norm(fixed),eps);
context.FixedPointCorrection = fixed-rawFixedResponse;
phiCritical = @(x) l6ns_phi(x,wc,context);
compensatedFixedResidual = norm(phiCritical(fixed)-fixed) / max(norm(fixed),eps);

% Verify that the nonlinear response map and saved analytic Jacobian use the
% same L6/IKp parameters and derivative path at the critical state.
rng(cfg.RandomSeed);
jvpError=nan(4,1);
jvpStep=2e-6*max(1,norm(fixed)/sqrt(numel(fixed)));
for directionIndex=1:numel(jvpError)
    direction=randn(size(fixed)); direction=direction/norm(direction);
    numeric=(phiCritical(fixed+jvpStep*direction)-phiCritical(fixed-jvpStep*direction))/(2*jvpStep);
    analytic=jCritical*direction;
    jvpError(directionIndex)=norm(numeric-analytic)/max(norm(analytic),eps);
end
jvpTable=table((1:numel(jvpError))',jvpError, ...
    'VariableNames',{'direction','relativeJvpError'});
writetable(jvpTable,fullfile(sectionDir,'analytic_numeric_jvp_check.tsv'), ...
    'FileType','text','Delimiter','\t');

stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
epsilon = stateScale * [1e-2 3e-3 1e-3 3e-4 1e-4];
quadratic = nan(numel(epsilon),1);
cubic = nan(numel(epsilon),1);
for ei = 1:numel(epsilon)
    h = epsilon(ei);
    phi0 = phiCritical(fixed);
    phiPlus = phiCritical(fixed+h*r);
    phiMinus = phiCritical(fixed-h*r);
    secondDirectional = (phiPlus-2*phi0+phiMinus)/(h^2);
    quadratic(ei) = 0.5*l'*secondDirectional;

    phiPlus2 = phiCritical(fixed+2*h*r);
    phiMinus2 = phiCritical(fixed-2*h*r);
    thirdDirectional = (phiPlus2-2*phiPlus+2*phiMinus-phiMinus2)/(2*h^3);
    cubic(ei) = (1/6)*l'*thirdDirectional;
end
normalFormTable = table(epsilon(:),quadratic,cubic, ...
    'VariableNames',{'epsilon','quadraticCoefficient','cubicCoefficient'});
writetable(normalFormTable,fullfile(sectionDir,'normal_form_convergence.tsv'), ...
    'FileType','text','Delimiter','\t');

validQuadratic = quadratic(isfinite(quadratic));
[b,quadraticRelativeSpread] = local_plateau_estimate(validQuadratic);
[cubicEstimate,cubicRelativeSpread] = local_plateau_estimate(cubic(isfinite(cubic)));
critical = struct('W',wc,'Alpha',alphaC,'Right',r,'Left',l, ...
    'Beta',beta,'Quadratic',b);

% Nonlinear perturbation validation in the original, fixed, and dominant
% compensated block-ablated systems.
[~,dominantIndex] = max(abs(section1.BlockSensitivity.sensitivityReal(2:end)));
dominantIndex = dominantIndex + 1;
dominantBlock = char(section1.BlockSensitivity.block(dominantIndex));
validationTable = local_nonlinear_validation(context,critical,dominantBlock,stateScale);
writetable(validationTable,fullfile(sectionDir,'nonlinear_perturbation_validation.tsv'), ...
    'FileType','text','Delimiter','\t');

% Continue a possible second branch through the persistent fixed point.
continuation = l6ns_continue_branch(context,critical,cfg);
continuation.Table.persistentMaxReal=nan(height(continuation.Table),1);
for rowIndex=1:height(continuation.Table)
    continuation.Table.persistentMaxReal(rowIndex)=l6ns_max_real( ...
        data.A+(1-continuation.Table.w(rowIndex))*data.B,cfg);
end
writetable(continuation.Table,fullfile(sectionDir,'second_branch_continuation.tsv'), ...
    'FileType','text','Delimiter','\t');

% Check the zero state, clamp margins, and active-set changes on the branch.
rawZeroResidual = norm(l6ns_phi(zeros(size(fixed)),wc,rawContext));
compensatedZeroResidual = norm(l6ns_phi(zeros(size(fixed)),wc,context));
[~,fixedL6Unclamped]=l6ns_l6_dynamic(fixed,context,[1 1]);
lowerMargin = min(fixedL6Unclamped-1);
upperMargin = min(40-fixedL6Unclamped);
minimumClampMargin = min(lowerMargin,upperMargin);
nearClampCount = sum(min(abs(fixedL6Unclamped-1),abs(40-fixedL6Unclamped))<1e-3);
fixedActive=fixedL6Unclamped>1 & fixedL6Unclamped<40;
clampRows=cell(height(continuation.Table),1);
for rowIndex=1:height(continuation.Table)
    state=continuation.States{rowIndex};
    [~,unclamped]=l6ns_l6_dynamic(state,context,[1 1]);
    active=unclamped>1 & unclamped<40;
    clampRows{rowIndex}={continuation.Table.amplitude(rowIndex),continuation.Table.w(rowIndex), ...
        min(abs(unclamped-1)),min(abs(40-unclamped)),sum(active~=fixedActive)};
end
clampTable=cell2table(vertcat(clampRows{:}),'VariableNames', ...
    {'amplitude','w','minimumDistanceToLowerClamp','minimumDistanceToUpperClamp','activeSetChanges'});
writetable(clampTable,fullfile(sectionDir,'continued_branch_clamp_check.tsv'), ...
    'FileType','text','Delimiter','\t');

simpleRealCrossing = abs(imag(lambdaCritical)) < 1e-7 && spectralSeparation > 1e-5 ...
    && secondLeadingReal < 1-1e-6;
transverseCrossing = abs(real(beta)) > 1e-5;
quadraticNonzero = abs(b) > 1e-5 && quadraticRelativeSpread < 0.35;
cubicNonzero = abs(cubicEstimate) > 1e-8 && cubicRelativeSpread < 0.35;
smoothAtClamp = minimumClampMargin > 1e-3;
continued = sum(continuation.Table.converged) >= 4;
bothSigns = any(continuation.Table.converged & continuation.Table.amplitude<0) ...
    && any(continuation.Table.converged & continuation.Table.amplitude>0);
validContinuation=continuation.Table.converged & isfinite(continuation.Table.maxRealLambda);
persistentStable=continuation.Table.persistentMaxReal<1;
otherStable=continuation.Table.maxRealLambda<1;
stabilityExchange=any(validContinuation & persistentStable & ~otherStable) && ...
    any(validContinuation & ~persistentStable & otherStable);
activeSetUnchanged=all(clampTable.activeSetChanges(continuation.Table.converged)==0);
jvpAligned=max(jvpError)<5e-3;
if ~jvpAligned
    classification = "analytic_nonlinear_derivative_mismatch";
elseif simpleRealCrossing && transverseCrossing && quadraticNonzero && smoothAtClamp && ...
        continued && bothSigns && stabilityExchange && activeSetUnchanged && jvpAligned
    classification = "transcritical_supported";
elseif simpleRealCrossing && transverseCrossing && quadraticNonzero && smoothAtClamp
    classification = "local_transcritical_coefficients_but_branch_not_fully_continued";
elseif simpleRealCrossing && transverseCrossing && ~quadraticNonzero && cubicNonzero && ...
        smoothAtClamp && continued && bothSigns && stabilityExchange && activeSetUnchanged
    classification = "pitchfork_like_supported";
elseif simpleRealCrossing && transverseCrossing && ~quadraticNonzero && cubicNonzero && smoothAtClamp
    classification = "pitchfork_or_symmetry_driven_candidate";
elseif ~smoothAtClamp
    classification = "possible_nonsmooth_or_border_collision";
elseif ~simpleRealCrossing
    classification = "not_a_simple_real_crossing";
elseif ~quadraticNonzero
    classification = "quadratic_small_pitchfork_or_degenerate_candidate";
else
    classification = "inconclusive";
end

summary = table(cfg.Angle,cfg.Contrast,wc,alphaC,real(lambdaCritical),imag(lambdaCritical), ...
    spectralSeparation,secondLeadingReal,real(beta),imag(beta),b,cubicEstimate, ...
    quadraticRelativeSpread,cubicRelativeSpread, ...
    rawFixedResidual,compensatedFixedResidual,fixedPointDifference, ...
    rawZeroResidual,compensatedZeroResidual,lowerMargin,upperMargin, ...
    minimumClampMargin,nearClampCount,simpleRealCrossing,transverseCrossing, ...
    quadraticNonzero,cubicNonzero,smoothAtClamp,continued,bothSigns,stabilityExchange, ...
    activeSetUnchanged,max(jvpError),jvpAligned,classification, ...
    'VariableNames',{'angle','contrast','wCritical','alphaCritical','lambdaReal','lambdaImag', ...
    'spectralSeparation','secondLeadingReal','betaReal','betaImag', ...
    'quadraticCoefficient','cubicCoefficient','quadraticRelativeSpread','cubicRelativeSpread', ...
    'rawFixedPointResidual','compensatedFixedPointResidual', ...
    'fixedPointEndpointDifference','rawZeroStateResidual','compensatedZeroStateResidual', ...
    'lowerClampMargin','upperClampMargin','minimumClampMargin','nearClampPixelCount', ...
    'simpleRealCrossing','transverseCrossing','quadraticNonzero','cubicNonzero','smoothAtClamp', ...
    'continuedBranch','continuedBothSigns','stabilityExchange','activeSetUnchanged', ...
    'maximumRelativeJvpError','analyticNumericJvpAligned','classification'});
writetable(summary,fullfile(sectionDir,'section3_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

local_plot_normal_form(normalFormTable,sectionDir);
local_plot_validation(validationTable,wc,sectionDir);
local_plot_continuation(continuation.Table,wc,sectionDir);
local_plot_branch_difference(continuation,data,sectionDir);

result = struct('Summary',summary,'NormalForm',normalFormTable, ...
    'Validation',validationTable,'Continuation',continuation, ...
    'ClampCheck',clampTable,'JvpCheck',jvpTable,'CriticalRight',r,'CriticalLeft',l);
save(fullfile(sectionDir,'section3_result.mat'),'result','-v7.3');
end

function tableOut = local_nonlinear_validation(context,critical,dominantBlock,stateScale)
weights = [critical.W+0.02,critical.W-0.02];
epsilonValues = stateScale*[1e-4 5e-4 1e-3];
signValues = [-1 1];
systems = {'full','fixed_L6',['remove_' dominantBlock]};
steps = 45;
rows = {};
for wi = 1:numel(weights)
    w = weights(wi);
    for systemIndex = 1:numel(systems)
        systemName = systems{systemIndex};
        [targetWeight,sourceScale] = local_system_control(w,systemName,dominantBlock);
        baselineTrajectory = cell(steps+1,1);
        baselineTrajectory{1} = context.FixedPoint;
        for baselineStep = 1:steps
            baselinePhi = l6ns_phi(baselineTrajectory{baselineStep}, ...
                targetWeight,context,sourceScale);
            baselineTrajectory{baselineStep+1} = baselineTrajectory{baselineStep} + ...
                context.RelaxationP*(baselinePhi-baselineTrajectory{baselineStep});
        end
        for ei = 1:numel(epsilonValues)
            for si = 1:numel(signValues)
                epsilon = signValues(si)*epsilonValues(ei);
                state = context.FixedPoint+epsilon*critical.Right;
                for step = 0:steps
                    baselineState = baselineTrajectory{step+1};
                    displacement = state-baselineState;
                    amplitude = critical.Left'*displacement;
                    rows(end+1,:) = {w,systemName,epsilon,step,norm(displacement), ...
                        norm(baselineState-context.FixedPoint), ...
                        real(amplitude),imag(amplitude),targetWeight(1),targetWeight(2), ...
                        targetWeight(3),sourceScale(1),sourceScale(2)}; %#ok<AGROW>
                    if step < steps
                        phi = l6ns_phi(state,targetWeight,context,sourceScale);
                        state = state+context.RelaxationP*(phi-state);
                    end
                end
            end
        end
    end
end
tableOut = cell2table(rows,'VariableNames',{'w','system','epsilon','step', ...
    'perturbationNorm','baselineDriftNorm','leftAmplitudeReal','leftAmplitudeImag', ...
    'wS','wC','wI','sourceScaleS','sourceScaleC'});
end

function [targetWeight,sourceScale] = local_system_control(w,systemName,dominantBlock)
targetWeight = [w w w];
sourceScale = [1 1];
if strcmp(systemName,'fixed_L6')
    targetWeight = [1 1 1];
elseif startsWith(systemName,'remove_target_')
    label = extractAfter(systemName,'remove_target_');
    index = find(strcmp({'S','C','I'},label),1);
    targetWeight(index) = 1;
elseif startsWith(systemName,'remove_source_')
    label = extractAfter(systemName,'remove_source_');
    index = find(strcmp({'S','C'},label),1);
    sourceScale(index) = 0;
elseif startsWith(systemName,'remove_')
    % The caller supplies names such as target_S or source_C.
    if startsWith(dominantBlock,'target_')
        label = extractAfter(dominantBlock,'target_');
        index = find(strcmp({'S','C','I'},label),1);
        targetWeight(index) = 1;
    elseif startsWith(dominantBlock,'source_')
        label = extractAfter(dominantBlock,'source_');
        index = find(strcmp({'S','C'},label),1);
        sourceScale(index) = 0;
    end
end
end

function [estimate,relativeSpread] = local_plateau_estimate(values)
values = values(:);
if isempty(values)
    estimate = NaN;
    relativeSpread = NaN;
    return
end
if numel(values) < 3
    estimate = median(values);
    relativeSpread = std(values)/max(abs(estimate),1e-12);
    return
end
score = nan(numel(values)-2,1);
windowMedian = nan(size(score));
for index = 1:numel(score)
    window = values(index:index+2);
    windowMedian(index) = median(window);
    score(index) = std(window)/max(abs(windowMedian(index)),1e-12);
end
[relativeSpread,index] = min(score);
estimate = windowMedian(index);
end

function local_plot_normal_form(t,outputDir)
fig = figure('Visible','off','Color','w','Position',[100 100 900 420]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; semilogx(t.epsilon,t.quadraticCoefficient,'-o'); grid on;
xlabel('\epsilon'); ylabel('b(\epsilon)'); title('quadratic coefficient convergence');
nexttile; semilogx(t.epsilon,t.cubicCoefficient,'-o'); grid on;
xlabel('\epsilon'); ylabel('c(\epsilon)'); title('cubic coefficient convergence');
l6ns_save_figure(fig,outputDir,'normal_form_coefficients'); close(fig);
end

function local_plot_validation(t,wc,outputDir)
fig = figure('Visible','off','Color','w','Position',[100 100 1250 520]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for wi = 1:2
    nexttile; hold on;
    wValues = unique(t.w);
    w = wValues(wi);
    subset = t(t.w==w & t.epsilon>0,:);
    systems = unique(subset.system,'stable');
    for si = 1:numel(systems)
        one = subset(strcmp(subset.system,systems{si}) & subset.epsilon==max(subset.epsilon),:);
        semilogy(one.step,one.perturbationNorm,'LineWidth',1.4);
    end
    grid on; xlabel('damped iteration step'); ylabel('||f_n-f_*||');
    title(sprintf('w=%.5f (%s threshold)',w,local_side(w,wc)));
    legend(systems,'Interpreter','none','Location','best');
end
l6ns_save_figure(fig,outputDir,'nonlinear_growth_decay'); close(fig);
end

function label = local_side(w,wc)
if w>wc; label='stable side'; else; label='unstable side'; end
end

function local_plot_continuation(t,wc,outputDir)
fig = figure('Visible','off','Color','w','Position',[100 100 1000 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; scatter(t.w,t.amplitude,55,double(t.converged),'filled'); hold on;
xline(wc,'k--'); grid on; xlabel('w'); ylabel('a=l^*(f-f_*)'); title('continued branch');
nexttile; scatter(t.w,t.maxRealLambda,55,double(t.converged),'filled'); hold on;
yline(1,'r--'); xline(wc,'k--'); grid on; xlabel('w'); ylabel('max Re \lambda'); title('second-branch stability');
l6ns_save_figure(fig,outputDir,'second_branch_continuation'); close(fig);
end

function local_plot_branch_difference(continuation,data,outputDir)
valid = find(continuation.Table.converged,1,'last');
if isempty(valid); return; end
state = continuation.States{valid};
if isempty(state); return; end
state = state-continuation.FixedPoint;
n = data.PopulationSize;
maps = {state(1:n),state(n+(1:n)),state(2*n+(1:n))};
populationLabels={'S','C','I'};
fig = figure('Visible','off','Color','w','Position',[100 100 1200 360]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i=1:3
    nexttile; imagesc(reshape(real(maps{i}),data.MapSize)); axis image off; colorbar;
    title(sprintf('%s branch difference',populationLabels{i}));
end
l6ns_save_figure(fig,outputDir,'second_branch_representative_state'); close(fig);
end
