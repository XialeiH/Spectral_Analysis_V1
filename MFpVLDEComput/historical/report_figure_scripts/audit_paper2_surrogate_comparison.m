function audit_paper2_surrogate_comparison()
% Reproduce the Paper2 pixelwise comparison on the newly simulated SNN data.

projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project/NYU-EffiModel';
artifactRoot = ['/Users/xialeihuang/Desktop/Neuroscience_Project/Spectral_Analysis/' ...
    'matlab-inserting_into_CG_model/report_figure_artifacts/' ...
    '1_Neural_Network_Surrogate'];
simulationRoot = fullfile(artifactRoot, 'snn_simulations');
dataRoot = fullfile(projectRoot, 'Data', 'Paper2_NetworkTuning', 'Fig1V4');
functionRoot = fullfile(dataRoot, 'Paper2PlotingData');
addpath(fullfile(projectRoot, 'Utils'), '-begin');

parameterData = load(fullfile(dataRoot, ...
    'AllMFPixPara_Paper2TuneFig1V4D2.mat'), 'CplxR');
angles = [0, 7.5, 15, 22.5];

networkE = [];
networkI = [];
paper2E = [];
paper2I = [];
angleDeg = [];
replicateIndex = [];
for ai = 1:numel(angles)
    responseFile = fullfile(functionRoot, sprintf( ...
        'Func16V4D2Ang%.1fSmall.mat', angles(ai)));
    responseData = load(responseFile, 'L4EmeshX', 'L4ImeshY', 'LDEFrfunc');
    for replicate = 1:2
        simulationFile = fullfile(simulationRoot, sprintf( ...
            'NWSimulationPix_%.1fdeg_rep%02d_h96_validation.mat', ...
            angles(ai), replicate));
        simulationData = load(simulationFile, 'NWSmlt');
        N = simulationData.NWSmlt;

        responseS = LDEIterFunc_Grating_16Func(responseData.L4EmeshX, ...
            responseData.L4ImeshY, responseData.LDEFrfunc.S, ...
            N.L4E, N.L4I, N.PixInptCtgrUse);
        responseC = LDEIterFunc_Grating_16Func(responseData.L4EmeshX, ...
            responseData.L4ImeshY, responseData.LDEFrfunc.C, ...
            N.L4E, N.L4I, N.PixInptCtgrUse);
        responseI = LDEIterFunc_Grating_16Func(responseData.L4EmeshX, ...
            responseData.L4ImeshY, responseData.LDEFrfunc.I, ...
            N.L4E, N.L4I, N.PixInptCtgrUse);

        networkEOne = N.FS * (1 - parameterData.CplxR) + ...
            N.FC * parameterData.CplxR;
        responseEOne = InhKill(responseS * (1 - parameterData.CplxR) + ...
            responseC * parameterData.CplxR, 50, [200 150], 200, 0);
        responseIOne = InhKill(responseI, 70, [100 95], 120, 0.9);

        networkE = [networkE; networkEOne]; %#ok<AGROW>
        networkI = [networkI; N.FI]; %#ok<AGROW>
        paper2E = [paper2E; responseEOne]; %#ok<AGROW>
        paper2I = [paper2I; responseIOne]; %#ok<AGROW>
        angleDeg = [angleDeg; repmat(angles(ai), numel(N.FI), 1)]; %#ok<AGROW>
        replicateIndex = [replicateIndex; repmat(replicate, numel(N.FI), 1)]; %#ok<AGROW>
    end
end

currentData = load(fullfile(artifactRoot, ...
    '1.2_1.3_real_paired_data_source.mat'), 'paired');
P = currentData.paired;

paper2 = summarizeComparison(networkE, paper2E, networkI, paper2I);
current = summarizeComparison(P.E_premodel, P.E_h96, ...
    P.I_premodel, P.I_h96);

fprintf('\nExact Paper2 library/postprocessing on new SNN simulations (7200 points):\n');
printSummary(paper2);
fprintf('\nCurrent corrected-h96 conversion (4680 points):\n');
printSummary(current);

byAngle = zeros(numel(angles), 8);
for ai = 1:numel(angles)
    mask = angleDeg == angles(ai);
    one = summarizeComparison(networkE(mask), paper2E(mask), ...
        networkI(mask), paper2I(mask));
    byAngle(ai,:) = [one.E.within33, one.I.within33, ...
        one.E.within25Ratio, one.I.within25Ratio, ...
        one.E.correlation, one.I.correlation, ...
        one.E.medianRelativeError, one.I.medianRelativeError];
end

paper2Paired = struct('networkE',networkE,'modelE',paper2E, ...
    'networkI',networkI,'modelI',paper2I,'angleDeg',angleDeg, ...
    'replicate',replicateIndex);
diagnostic = struct('paper2', paper2, 'currentH96', current, ...
    'anglesDeg', angles, 'paper2ByAngle', byAngle, ...
    'paper2ByAngleColumns', {{'E_within33','I_within33', ...
    'E_within_ratio_0p8_1p25','I_within_ratio_0p8_1p25', ...
    'E_correlation','I_correlation','E_median_relative_error', ...
    'I_median_relative_error'}});
save(fullfile(artifactRoot, '1.2_1.3_paper2_reproduction_audit.mat'), ...
    'diagnostic', 'paper2Paired', '-v7.3');
plotPaper2Reproduction(artifactRoot, paper2Paired, angles);

reportFile = fullfile(artifactRoot, ...
    '1.2_1.3_paper2_reproduction_audit.txt');
fid = fopen(reportFile, 'w');
assert(fid >= 0, 'Cannot create %s.', reportFile);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
writeSummary(fid, 'Exact Paper2 library/postprocessing', paper2);
writeSummary(fid, 'Current corrected-h96 conversion', current);
fprintf(fid, '\nPaper2 reproduction by angle:\n');
fprintf(fid, ['angle  E_33  I_33  E_ratio_0.8_1.25  I_ratio_0.8_1.25  ' ...
    'E_corr  I_corr  E_median_relerr  I_median_relerr\n']);
for ai = 1:numel(angles)
    fprintf(fid, '%.1f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f  %.8f\n', ...
        angles(ai), byAngle(ai,:));
end

function plotPaper2Reproduction(outputRoot,P,angles)
% Match the original Paper2 display density: 25% of one replicate.
keep = false(size(P.networkE));
for ai = 1:numel(angles)
    candidates = find(P.angleDeg == angles(ai) & P.replicate == 1);
    selected = unique(round(linspace(1,numel(candidates),225)));
    keep(candidates(selected)) = true;
end

fig = figure('Visible','off','Color','w','Units','inches', ...
    'Position',[1 1 13.2 6.1]);
t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(t);
drawPanel(ax,P.networkE(keep),P.modelE(keep),[5 40],'E',[1 0.125 0.125]);
xticks(ax,[10 20 30 40]); yticks(ax,5:5:40);
ax = nexttile(t);
drawPanel(ax,P.networkI(keep),P.modelI(keep),[25 110],'I',[0.21 0.36 1]);
xticks(ax,[40 60 80 100]); yticks(ax,30:10:110);

base = fullfile(outputRoot,'diagnostic_Paper2_exact_reproduction');
savefig(fig,[base '.fig']);
exportgraphics(fig,[base '.png'],'Resolution',600);
exportgraphics(fig,[base '.pdf'],'ContentType','vector');
close(fig);
end

function drawPanel(ax,x,y,limits,population,color)
hold(ax,'on'); lineX=linspace(limits(1),limits(2),500);
plot(ax,lineX,lineX,'k--','LineWidth',1.1);
plot(ax,lineX,(4/3)*lineX,'k--','LineWidth',1.5);
plot(ax,lineX,lineX/(4/3),'k--','LineWidth',1.5);
scatter(ax,x,y,1,'o','filled','MarkerFaceColor',color,'MarkerEdgeColor','none');
xlim(ax,limits); ylim(ax,limits); axis(ax,'square');
title(ax,sprintf('Pixelwise %s firing rates',population), ...
    'FontSize',20,'FontWeight','bold');
xlabel(ax,'Premodel Fr (Hz)','FontSize',18);
if population=='E'; ylabel(ax,'Paper2 CG Model Fr (Hz)','FontSize',18); end
set(ax,'FontName','Arial','FontSize',15,'LineWidth',1.2,'Box','off', ...
    'TickDir','in','Layer','top');
end
fprintf('Saved audit to %s\n', reportFile);
end

function result = summarizeComparison(networkE, modelE, networkI, modelI)
result = struct('E', populationSummary(networkE, modelE), ...
    'I', populationSummary(networkI, modelI));
end

function result = populationSummary(network, model)
valid = isfinite(network) & isfinite(model) & network > 0;
network = network(valid);
model = model(valid);
ratio = model ./ network;
relativeError = abs(model - network) ./ network;
result = struct('count', numel(network), ...
    'within33', mean(relativeError <= 0.33), ...
    'within20', mean(relativeError <= 0.20), ...
    'within25Ratio', mean(ratio >= 0.8 & ratio <= 1.25), ...
    'medianRelativeError', median(relativeError), ...
    'meanRelativeError', mean(relativeError), ...
    'correlation', corr(network, model), ...
    'networkMean', mean(network), 'modelMean', mean(model));
end

function printSummary(result)
fprintf(['E: n=%d, within33=%.2f%%, within20=%.2f%%, ratio[0.8,1.25]=%.2f%%, ' ...
    'median relerr=%.4f, corr=%.4f, means=%.3f/%.3f\n'], ...
    result.E.count, 100*result.E.within33, 100*result.E.within20, ...
    100*result.E.within25Ratio, result.E.medianRelativeError, ...
    result.E.correlation, result.E.networkMean, result.E.modelMean);
fprintf(['I: n=%d, within33=%.2f%%, within20=%.2f%%, ratio[0.8,1.25]=%.2f%%, ' ...
    'median relerr=%.4f, corr=%.4f, means=%.3f/%.3f\n'], ...
    result.I.count, 100*result.I.within33, 100*result.I.within20, ...
    100*result.I.within25Ratio, result.I.medianRelativeError, ...
    result.I.correlation, result.I.networkMean, result.I.modelMean);
end

function writeSummary(fid, label, result)
fprintf(fid, '\n%s\n', label);
fprintf(fid, ['E count=%d within33=%.8f within20=%.8f ratio_0p8_1p25=%.8f ' ...
    'median_relative_error=%.8f mean_relative_error=%.8f correlation=%.8f ' ...
    'network_mean=%.8f model_mean=%.8f\n'], result.E.count, ...
    result.E.within33, result.E.within20, result.E.within25Ratio, ...
    result.E.medianRelativeError, result.E.meanRelativeError, ...
    result.E.correlation, result.E.networkMean, result.E.modelMean);
fprintf(fid, ['I count=%d within33=%.8f within20=%.8f ratio_0p8_1p25=%.8f ' ...
    'median_relative_error=%.8f mean_relative_error=%.8f correlation=%.8f ' ...
    'network_mean=%.8f model_mean=%.8f\n'], result.I.count, ...
    result.I.within33, result.I.within20, result.I.within25Ratio, ...
    result.I.medianRelativeError, result.I.meanRelativeError, ...
    result.I.correlation, result.I.networkMean, result.I.modelMean);
end
