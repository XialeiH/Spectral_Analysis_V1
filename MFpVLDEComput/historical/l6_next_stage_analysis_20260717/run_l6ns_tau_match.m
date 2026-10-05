% Match the nonlinear ODE time scale to the canonical relaxed iteration.

cfg = l6ns_config();
data = l6ns_load_endpoints(cfg);
canonicalScript = getenv('L6NS_CANONICAL_MODEL_SCRIPT');
if isempty(canonicalScript)
    error('L6NS_CANONICAL_MODEL_SCRIPT must name the canonical h96 setup script.');
end
if isempty(getenv('L6NS_CANONICAL_CODE_DIR'))
    error('L6NS_CANONICAL_CODE_DIR must name the canonical Torch geometry code directory.');
end
setenv('GEOM_ANGLE_COND',sprintf('%.12g',cfg.Angle));
setenv('GEOM_CONTR_COND',sprintf('%.12g',cfg.Contrast));
setenv('L6_EQU_BLEND_WEIGHT','0');
setenv('L6NS_SETUP_ONLY','1');
run(canonicalScript);
run(fullfile(cfg.CodeRoot,'l6ns_context_from_canonical_workspace.m'));
tauMatch = l6ns_match_tau_iteration_ode(cfg,data,ModelContext);
