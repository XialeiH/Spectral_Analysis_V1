% Slurm entry point for near-silent inhibition branch landmarks.

setupFile = getenv('L6NS_GLOBAL_BIFURCATION_SETUP');
branchFile = getenv('L6NS_INHIBITION_NEARSILENT_BRANCH');
outputDir = getenv('L6NS_GLOBAL_BIFURCATION_EXTRA_ROOT');
if isempty(setupFile) || isempty(branchFile) || isempty(outputDir)
    error('Setup, near-silent branch, and output roots are required.');
end
setupLoaded = load(setupFile,'setup');
branchLoaded = load(branchFile,'result');
tableOut = l6ns_inhibition_nearsilent_landmarks(setupLoaded.setup, ...
    branchLoaded.result,outputDir); %#ok<NASGU>
