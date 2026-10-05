% Run Section 5 after Sections 1 and 2 have completed.
cfg=l6ns_config();
data=l6ns_load_endpoints(cfg);
loaded1=load(fullfile(cfg.OutputRoot,'section1_l6_causality','section1_result.mat'),'result');
loaded2=load(fullfile(cfg.OutputRoot,'section2_shared_modes','section2_result.mat'),'result');
section5=l6ns_section5_reduced_subspaces(cfg,data,loaded1.result,loaded2.result);
