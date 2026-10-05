% Run Section 3 after executing the canonical one-condition h96 setup.

cfg = l6ns_config();
data = l6ns_load_endpoints(cfg);
section1File = fullfile(cfg.OutputRoot,'section1_l6_causality','section1_result.mat');
assert(isfile(section1File),'Run Section 1 first: %s is missing.',section1File);
loadedSection1 = load(section1File,'result');
section1 = loadedSection1.result;

canonicalScript = getenv('L6NS_CANONICAL_MODEL_SCRIPT');
if isempty(canonicalScript)
    error('L6NS_CANONICAL_MODEL_SCRIPT must name the canonical h96 setup script.');
end
setenv('GEOM_ANGLE_COND',sprintf('%.12g',cfg.Angle));
setenv('GEOM_CONTR_COND',sprintf('%.12g',cfg.Contrast));
setenv('L6_EQU_BLEND_WEIGHT','0');
setenv('L6NS_SETUP_ONLY','1');
if isempty(getenv('L6NS_CANONICAL_CODE_DIR'))
    error('L6NS_CANONICAL_CODE_DIR must name the canonical Torch geometry code directory.');
end
run(canonicalScript);
run(fullfile(cfg.CodeRoot,'l6ns_context_from_canonical_workspace.m'));
section3 = l6ns_section3_bifurcation(cfg,data,section1,ModelContext);
