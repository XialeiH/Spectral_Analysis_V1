% Run Section 2 after Section 1 has completed.
cfg=l6ns_config();
data=l6ns_load_endpoints(cfg);
loaded=load(fullfile(cfg.OutputRoot,'section1_l6_causality','section1_result.mat'),'result');
section2=l6ns_section2_shared_modes(cfg,data,loaded.result);
