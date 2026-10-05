% Slurm entry point for the inhibition near-silent fixed-point branch.

setupFile = getenv('L6NS_GLOBAL_BIFURCATION_SETUP');
multistartFile = getenv('L6NS_INHIBITION_W1_MULTISTART');
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_EXTRA_ROOT');
if isempty(setupFile) || isempty(multistartFile) || isempty(outputDir)
    error('Setup, inhibition w=1 multistart, and extra output roots are required.');
end
setupLoaded = load(setupFile,'setup');
multistartLoaded = load(multistartFile,'result');
result = l6ns_inhibition_nearsilent_branch(setupLoaded.setup, ...
    multistartLoaded.result,outputDir); %#ok<NASGU>
