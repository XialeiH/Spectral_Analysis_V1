% Prepare the reusable canonical context for global bifurcation experiments.

cfg = l6ns_config();
data = l6ns_load_endpoints(cfg); %#ok<NASGU>
canonicalScript = getenv('L6NS_CANONICAL_MODEL_SCRIPT');
if isempty(canonicalScript)
    error('L6NS_CANONICAL_MODEL_SCRIPT is required.');
end
twoPathwayFile = getenv('L6NS_TWO_PATHWAY_RESULT');
if isempty(twoPathwayFile)
    error('L6NS_TWO_PATHWAY_RESULT is required.');
end
setenv('GEOM_ANGLE_COND',sprintf('%.12g',cfg.Angle));
setenv('GEOM_CONTR_COND',sprintf('%.12g',cfg.Contrast));
setenv('L6_EQU_BLEND_WEIGHT','0');
setenv('L6NS_SETUP_ONLY','1');
run(canonicalScript);
run(fullfile(cfg.CodeRoot,'l6ns_context_from_canonical_workspace.m'));
loaded = load(twoPathwayFile,'result');
setup = l6ns_global_bifurcation_setup(cfg,ModelContext,loaded.result.Pathway); %#ok<NASGU>
