% Run paired nonlinear bifurcation tests for L6 and inhibition.

cfg = l6ns_config();
data = l6ns_load_endpoints(cfg);
canonicalScript = getenv('L6NS_CANONICAL_MODEL_SCRIPT');
if isempty(canonicalScript)
    error('L6NS_CANONICAL_MODEL_SCRIPT must name the canonical h96 setup script.');
end
setenv('GEOM_ANGLE_COND',sprintf('%.12g',cfg.Angle));
setenv('GEOM_CONTR_COND',sprintf('%.12g',cfg.Contrast));
setenv('L6_EQU_BLEND_WEIGHT','0');
setenv('L6NS_SETUP_ONLY','1');
run(canonicalScript);
run(fullfile(cfg.CodeRoot,'l6ns_context_from_canonical_workspace.m'));

l6ResultFile = getenv('L6NS_L6_BIFURCATION_RESULT');
if isempty(l6ResultFile)
    error('L6NS_L6_BIFURCATION_RESULT must name section3_result.mat.');
end
twoPathwayFile = getenv('L6NS_TWO_PATHWAY_RESULT');
if isempty(twoPathwayFile)
    error('L6NS_TWO_PATHWAY_RESULT must name two_pathway_result.mat.');
end
l6Loaded = load(l6ResultFile,'result');
twoLoaded = load(twoPathwayFile,'result');
bifurcation = l6ns_compare_bifurcations( ...
    cfg,data,ModelContext,l6Loaded.result,twoLoaded.result);
