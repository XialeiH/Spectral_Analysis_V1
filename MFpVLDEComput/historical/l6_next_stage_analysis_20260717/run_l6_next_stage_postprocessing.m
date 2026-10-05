% Run Sections 1, 2, 4, and 5 from the saved endpoint Jacobians.

cfg=l6ns_config();
data=l6ns_load_endpoints(cfg);
section1=l6ns_section1_l6_causality(cfg,data);
section2=l6ns_section2_shared_modes(cfg,data,section1);
section4=l6ns_section4_component_modes(cfg,data,section1);
section5=l6ns_section5_reduced_subspaces(cfg,data,section1,section2);
save(fullfile(cfg.OutputRoot,'postprocessing_manifest.mat'),'cfg','-v7.3');
