function panel = l6ns_sparse_operating_panel(cfg,contrasts)
% Replicate pathway sensitivities across the available endpoint conditions.

outputDir = fullfile(cfg.OutputRoot,'sparse_operating_panel');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
rows = cell(numel(contrasts),1);
for index = 1:numel(contrasts)
    conditionCfg = cfg;
    conditionCfg.Contrast = contrasts(index);
    data = l6ns_load_endpoints(conditionCfg);
    pathway = l6ns_two_pathway_data(data);
    baseline = l6ns_pathway_jacobian(pathway,1,1);
    frozenL6 = l6ns_pathway_jacobian(pathway,0,1);
    frozenI = l6ns_pathway_jacobian(pathway,1,0);
    modes = l6ns_eigenpairs(baseline,6,'largestreal',conditionCfg);
    r = modes.Right(:,1); l = modes.Left(:,1);
    lambda = modes.Lambda(1);
    s6 = l'*pathway.J6*r;
    sI = l'*pathway.JI*r;
    alphaFrozenL6 = l6ns_max_real(frozenL6,conditionCfg);
    alphaFrozenI = l6ns_max_real(frozenI,conditionCfg);
    generator = (baseline-speye(data.Dimension))/cfg.TauMs;
    numerical = real(eigs((generator+generator')/2,1,'largestreal'));
    rows{index} = {data.Angle,data.Contrast,real(lambda),imag(lambda), ...
        (real(lambda)-1)/cfg.TauMs,cfg.TauMs/max(1-real(lambda),eps), ...
        real(s6),imag(s6),real(sI),imag(sI),alphaFrozenL6,alphaFrozenI, ...
        modes.ConditionNumber(1),1/modes.ConditionNumber(1),numerical, ...
        pathway.J6ISourceFraction,pathway.ReconstructionError};
end
summary = cell2table(vertcat(rows{:}),'VariableNames', ...
    {'angle','contrast','baselineLambdaReal','baselineLambdaImag', ...
    'baselineAlphaPhysicalPerMs','baselineRecoveryMs','s6Real','s6Imag', ...
    'sIReal','sIImag','frozenL6LambdaReal','frozenILambdaReal', ...
    'conditionNumber','leftRightAlignment','numericalAbscissaPerMs', ...
    'J6ISourceFraction','pathwayReconstructionError'});
writetable(summary,fullfile(outputDir,'sparse_operating_panel.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[100 100 1200 430]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; plot(summary.contrast,summary.baselineLambdaReal,'-o','LineWidth',1.4); hold on;
plot(summary.contrast,summary.frozenL6LambdaReal,'--o','LineWidth',1.4);
plot(summary.contrast,summary.frozenILambdaReal,':o','LineWidth',1.4); yline(1,'r--');
grid on; xlabel('contrast'); ylabel('max Re \lambda'); legend({'baseline','freeze L6','freeze I'});
nexttile; plot(summary.contrast,summary.s6Real,'-o','LineWidth',1.4); hold on;
plot(summary.contrast,summary.sIReal,'-o','LineWidth',1.4); grid on;
xlabel('contrast'); ylabel('pathway sensitivity'); legend({'Re(y^*J_6x)','Re(y^*J_Ix)'});
nexttile; semilogy(summary.contrast,summary.conditionNumber,'-o','LineWidth',1.4); grid on;
xlabel('contrast'); ylabel('critical condition number');
l6ns_save_figure(fig,outputDir,'sparse_operating_point_replication'); close(fig);
panel = struct('Summary',summary);
save(fullfile(outputDir,'sparse_operating_panel_result.mat'),'panel','-v7.3');
end
