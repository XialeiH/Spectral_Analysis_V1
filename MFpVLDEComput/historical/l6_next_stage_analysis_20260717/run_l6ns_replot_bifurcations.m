% Regenerate paired bifurcation figures from an existing result file.

resultFile = getenv('L6NS_BIFURCATION_RESULT');
if isempty(resultFile)
    error('L6NS_BIFURCATION_RESULT must name l6_inhibition_bifurcation_result.mat.');
end
outputDir = fileparts(resultFile);
loaded = load(resultFile,'result');
l6ns_plot_bifurcations(loaded.result,outputDir);
