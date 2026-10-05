function run_figure1g_densemul_smoke(runRoot)
% Verify dense-block multiplication against the sparse reference at 4x4.
addpath(fullfile(runRoot,'code'),'-begin');
run_figure1g_field_trial_densemul(runRoot,4,0);
end
