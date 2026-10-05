% Fit the current h96 ODE tau to the precomputed 0-degree LIF E peak.

cfg = l6ns_config();
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
fprintf('Resolved LocalResponse: %s\n',which('LocalResponse_6D_MLP_prefAngle'));
fprintf('Resolved L6Convert: %s\n',which('L6Convert'));
run(fullfile(cfg.CodeRoot,'l6ns_context_from_canonical_workspace.m'));
tauFit = l6ns_fit_tau_lif_peak(cfg,ModelContext);
