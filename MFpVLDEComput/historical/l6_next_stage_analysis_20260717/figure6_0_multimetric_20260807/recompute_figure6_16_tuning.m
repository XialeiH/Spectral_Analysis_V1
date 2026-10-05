function recompute_figure6_16_tuning()
% Recompute only the baseline and eight marked contrast-100 tuning curves.

datasetFile = getenv('FIG616_DATASET_FILE');
outputRoot = getenv('FIG616_CORRECTED_ROOT');
selectionIndex = str2double(getenv('FIG616_SELECTION_INDEX'));
if isempty(datasetFile) || isempty(outputRoot) || ...
        ~ismember(selectionIndex,0:8)
    error('Figure616:RecomputeEnvironment', ...
        ['FIG616_DATASET_FILE, FIG616_CORRECTED_ROOT, and ' ...
        'FIG616_SELECTION_INDEX=0:8 are required.']);
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

loaded = load(datasetFile,'setup');
setup = loaded.setup;
requestedBeta6 = [0.056 0.024 0.240 0.368 0.320, ...
                  0.168 0.168 0.176];
requestedBetaI = [0.016 0.104 0.096 0.112 0.352, ...
                  0.136 0.112 0.064];

if selectionIndex==0
    beta6 = 0;
    betaI = 0;
    canonicalStates = setup.BaselineCanonicalStates;
    canonicalResiduals = zeros(size(setup.CanonicalAngles));
else
    beta6 = requestedBeta6(selectionIndex);
    betaI = requestedBetaI(selectionIndex);
    gain6 = 1+beta6;
    gainI = 1+betaI;
    context = setup.Context;
    canonicalStates = zeros(size(setup.BaselineCanonicalStates));
    canonicalResiduals = nan(size(setup.CanonicalAngles));
    for angleIndex = 1:numel(setup.CanonicalAngles)
        angleContext = context;
        angleContext.OrientationUse = setup.CanonicalAngles(angleIndex)* ...
            ones(size(context.OrientationUse));
        phi = @(state)l6ns_phi(state,1-gain6,angleContext,[1 1],gainI,'true');
        fixed = real_tuning_fixed_point_tolerance( ...
            phi,setup.BaselineCanonicalStates(:,angleIndex), ...
            angleContext.RelaxationP,1e-4);
        canonicalResiduals(angleIndex) = fixed.Residual;
        if ~fixed.Converged
            error('Figure616:SelectedTuningFixedPoint', ...
                ['Selection %d, beta6 %.3f, betaI %.3f, angle %.2f ' ...
                'failed with residual %.6g.'],selectionIndex,beta6,betaI, ...
                setup.CanonicalAngles(angleIndex),fixed.Residual);
        end
        canonicalStates(:,angleIndex) = fixed.State;
    end
end

tuningCurves = figure6_0_tuning_curves( ...
    canonicalStates,setup.CanonicalAngles,setup.CWeight);
outputFile = fullfile(outputRoot,sprintf( ...
    'corrected_figure6_16_tuning_%d.mat',selectionIndex));
save(outputFile,'selectionIndex','beta6','betaI','tuningCurves', ...
    'canonicalResiduals','canonicalStates','-v7.3');
fprintf(['Corrected Figure 6.16 selection %d: beta6 %.3f betaI %.3f, ' ...
    'maximum residual %.6g -> %s\n'],selectionIndex,beta6,betaI, ...
    max(canonicalResiduals),outputFile);
end
