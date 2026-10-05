% Run Section 1 from saved w=0 and w=1 endpoint Jacobians.
cfg=l6ns_config();
data=l6ns_load_endpoints(cfg);
section1=l6ns_section1_l6_causality(cfg,data);
