% Run Section 4 after Section 1 has completed.
cfg=l6ns_config();
data=l6ns_load_endpoints(cfg);
loaded=load(fullfile(cfg.OutputRoot,'section1_l6_causality','section1_result.mat'),'result');
section4=l6ns_section4_component_modes(cfg,data,loaded.result);
