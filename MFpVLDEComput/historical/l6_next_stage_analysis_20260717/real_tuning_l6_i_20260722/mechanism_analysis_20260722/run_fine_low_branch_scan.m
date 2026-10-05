function run_fine_low_branch_scan()
% Resolve narrow stability changes on the raw-NN low-rate branches.

paths = mechanism_initialize();
loadedSetup = load(paths.SetupFile,'setup');
context = loadedSetup.setup.Context;
loaded = load(fullfile(paths.OutputRoot,'branch_continuation.mat'));

l6Coarse = 0:0.01:0.35;
l6Start = loaded.l6ForwardStates(:,find(abs(l6Coarse-0.16)<1e-12,1));
[l6Table,l6States] = local_scan("L6",0.16:0.001:0.18,l6Start,context);

iCoarse = 0:-0.005:-0.15;
iStart = loaded.iForwardStates(:,find(abs(iCoarse+0.07)<1e-12,1));
[iTable,iStates] = local_scan("Inhibition",-0.07:-0.001:-0.10,iStart,context);

fineScanTable = [l6Table;iTable];
writetable(fineScanTable,fullfile(paths.OutputRoot,'fine_low_branch_stability.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(paths.OutputRoot,'fine_low_branch_stability.mat'), ...
    'fineScanTable','l6States','iStates','-v7.3');

figure('Color','w','Position',[100 100 1180 460]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for index = 1:2
    nexttile;
    if index==1; rows=l6Table; else; rows=iTable; end
    plot(rows.beta,rows.maxRealLambda,'-o','LineWidth',1.5,'MarkerSize',4); hold on
    yline(1,'k--'); grid on
    xlabel('\beta'); ylabel('max Re(\lambda) of D\Phi');
    title(sprintf('%s low-rate branch',rows.pathway(1)));
end
title(layout,'Fine stability scan entirely inside the trained raw h96 NN domain');
exportgraphics(gcf,fullfile(paths.FigureRoot, ...
    '05_fine_low_branch_stability_scan.pdf'),'ContentType','vector');
end

function [rows,states] = local_scan(pathway,betaGrid,state,context)
cells = cell(numel(betaGrid),1);
states = nan(numel(state),numel(betaGrid));
for index = 1:numel(betaGrid)
    beta = betaGrid(index);
    gain6 = 1;
    gainI = 1;
    if pathway=="L6"; gain6=1+beta; else; gainI=1+beta; end
    phi = @(x)mechanism_phi_variant(x,context,gain6,gainI,'extended');
    fixedResult = real_tuning_fixed_point(phi,state,context.RelaxationP);
    state = fixedResult.State(:);
    states(:,index) = state;
    J = real_tuning_true_jacobian(state,context,gain6,gainI);
    maxReal = mechanism_max_real(J);
    stats = mechanism_domain_stats(state,context,gain6,gainI, ...
        pathway,beta,'fine low branch');
    cells{index} = table(pathway,beta,maxReal,fixedResult.Converged, ...
        fixedResult.Iterations,fixedResult.Residual,max(state), ...
        stats.fractionOutsideTrainingBox, ...
        'VariableNames',{'pathway','beta','maxRealLambda','converged', ...
        'iterations','residual','maximumRateHz','fractionOutsideTrainingBox'});
    fprintf('%s beta=%+.3f: maxRe %.6f, max rate %.3f, iter %d\n', ...
        pathway,beta,maxReal,max(state),fixedResult.Iterations);
end
rows = vertcat(cells{:});
end
