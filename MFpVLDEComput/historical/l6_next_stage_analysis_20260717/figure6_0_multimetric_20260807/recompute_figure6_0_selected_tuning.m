function recompute_figure6_0_selected_tuning()
% Recompute the three selected contrast-100 curves with the canonical remap.

runRoot = getenv('FIG60_RUN_ROOT');
selectionIndex = str2double(getenv('FIG60_SELECTED_INDEX'));
if isempty(runRoot) || ~ismember(selectionIndex,1:3)
    error('Figure60:SelectedTuningEnvironment', ...
        'FIG60_RUN_ROOT and FIG60_SELECTED_INDEX=1:3 are required.');
end

loaded = load(fullfile(runRoot,'setup', ...
    'figure6_0_baseline_contrast100.mat'),'setup');
setup = loaded.setup;
selectedBeta6 = [0.040 0.120 0.256];
selectedBetaI = [0.024 0.064 0.144];
beta6 = selectedBeta6(selectionIndex);
betaI = selectedBetaI(selectionIndex);
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
        error('Figure60:SelectedTuningFixedPoint', ...
            'Selection %d angle %.2f failed with residual %.6g.', ...
            selectionIndex,setup.CanonicalAngles(angleIndex),fixed.Residual);
    end
    canonicalStates(:,angleIndex) = fixed.State;
end

baselineTuning = figure6_0_tuning_curves( ...
    setup.BaselineCanonicalStates,setup.CanonicalAngles,setup.CWeight);
selectedTuning = figure6_0_tuning_curves( ...
    canonicalStates,setup.CanonicalAngles,setup.CWeight);
outputFile = fullfile(runRoot,'results',sprintf( ...
    'corrected_selected_tuning_%d.mat',selectionIndex));
save(outputFile,'selectionIndex','beta6','betaI','baselineTuning', ...
    'selectedTuning','canonicalResiduals','-v7.3');
fprintf(['Corrected selection %d: beta6 %.3f betaI %.3f, ' ...
    'maximum canonical residual %.6g -> %s\n'],selectionIndex,beta6,betaI, ...
    max(canonicalResiduals),outputFile);
end
