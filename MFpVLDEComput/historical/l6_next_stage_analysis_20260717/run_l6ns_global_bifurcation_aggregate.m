% Aggregate all wide-continuation and root-census experiments.

setupFile = getenv('L6NS_GLOBAL_BIFURCATION_SETUP');
branchRoot = getenv('L6NS_GLOBAL_BIFURCATION_BRANCH_ROOT');
multistartRoot = getenv('L6NS_GLOBAL_BIFURCATION_MULTISTART_ROOT');
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_FINAL_ROOT');
if isempty(setupFile) || isempty(branchRoot) || isempty(multistartRoot) || isempty(outputDir)
    error('Setup, branch, multistart, and final output roots are required.');
end
loaded = load(setupFile,'setup');
result = l6ns_global_bifurcation_aggregate(loaded.setup,branchRoot, ...
    multistartRoot,outputDir); %#ok<NASGU>
