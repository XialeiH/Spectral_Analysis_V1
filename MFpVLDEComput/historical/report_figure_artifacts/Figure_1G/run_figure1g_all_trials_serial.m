function run_figure1g_all_trials_serial(runRoot)
% Run every paired field-size benchmark without cross-job memory contention.
fields = [4 6 8 10 20 30 40];
for repeatId = 1:10
    for fieldHC = fields
        run_figure1g_field_trial(runRoot,fieldHC,repeatId);
    end
end
end
