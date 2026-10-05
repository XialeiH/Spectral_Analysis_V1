function aggregate_followup_results()
% Aggregate all follow-up tasks and create branch-resolved mechanism figures.

paths=followup_initialize();
branch=local_read(fullfile(paths.FollowupOutput,'branch_tasks','*.tsv'));
sensitivity=local_read(fullfile(paths.FollowupOutput,'sensitivity_tasks', ...
    'sensitivity_*.tsv'));
slopes=local_read(fullfile(paths.FollowupOutput,'sensitivity_tasks','slopes_*.tsv'));
interpolation=local_read(fullfile(paths.FollowupOutput,'interpolation_tasks', ...
    'interpolation_0*.tsv'));
factorial=local_read(fullfile(paths.FollowupOutput,'interpolation_tasks','factorial_*.tsv'));
shapley=local_read(fullfile(paths.FollowupOutput,'interpolation_tasks','shapley_*.tsv'));
meanPattern=local_read(fullfile(paths.FollowupOutput,'interpolation_tasks','mean_pattern_*.tsv'));
threshold=local_read(fullfile(paths.FollowupOutput,'threshold_tasks','threshold_*.tsv'));
odeTrajectory=readtable(fullfile(paths.FollowupOutput,'ode_hidden_mode', ...
    'ode_hidden_mode_trajectory.tsv'),'FileType','text','Delimiter','\t', ...
    'TextType','string');
odeSummary=readtable(fullfile(paths.FollowupOutput,'ode_hidden_mode', ...
    'ode_hidden_mode_summary.tsv'),'FileType','text','Delimiter','\t', ...
    'TextType','string');

writetable(sortrows(branch,{'pathway','beta'}), ...
    fullfile(paths.FollowupOutput,'branch_resolved_stability.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(sortrows(sensitivity,{'pathway','beta'}), ...
    fullfile(paths.FollowupOutput,'eigenvalue_sensitivity.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(slopes,fullfile(paths.FollowupOutput,'local_slope_distributions.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(interpolation,fullfile(paths.FollowupOutput,'state_interpolation.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(factorial,fullfile(paths.FollowupOutput,'population_factorial.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(shapley,fullfile(paths.FollowupOutput,'population_shapley.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(meanPattern,fullfile(paths.FollowupOutput,'mean_pattern_split.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(threshold,fullfile(paths.FollowupOutput,'extension_threshold_control.tsv'), ...
    'FileType','text','Delimiter','\t');

local_plot_branches(branch,paths.FollowupFigures);
local_plot_sensitivity(sensitivity,paths.FollowupFigures);
local_plot_interpolation(interpolation,shapley,meanPattern,paths.FollowupFigures);
local_plot_threshold(threshold,paths.FollowupFigures);
local_plot_ode(odeTrajectory,paths.FollowupFigures);
save(fullfile(paths.FollowupOutput,'followup_summary.mat'),'branch','sensitivity', ...
    'slopes','interpolation','factorial','shapley','meanPattern','threshold', ...
    'odeTrajectory','odeSummary','-v7.3');
fprintf('Aggregated follow-up investigation under %s.\n',paths.FollowupOutput);
end

function data=local_read(pattern)
files=dir(pattern);
if isempty(files); error('Followup:MissingFiles','No files match %s.',pattern); end
cells=cell(numel(files),1);
for index=1:numel(files)
    cells{index}=readtable(fullfile(files(index).folder,files(index).name), ...
        'FileType','text','Delimiter','\t','TextType','string');
end
data=vertcat(cells{:});
end

function local_plot_branches(data,figureRoot)
fig=figure('Visible','off','Color','w','Position',[100 100 1380 580]);
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
for p=1:2
    pathway=["L6","Inhibition"]; rows=data(data.pathway==pathway(p),:);
    rows=sortrows(rows,'beta'); ax=nexttile(layout); hold(ax,'on');
    plot(ax,rows.beta,rows.maxRealFpp,'-o','LineWidth',1.4,'DisplayName','FPP');
    plot(ax,rows.beta,rows.maxRealDirectBaseline,'-s','LineWidth',1.4, ...
        'DisplayName','true input at baseline');
    plot(ax,rows.beta,rows.maxRealLowBranch,'-^','LineWidth',1.8, ...
        'DisplayName','low moved branch');
    plot(ax,rows.beta,rows.maxRealHighBranch,'-d','LineWidth',1.8, ...
        'DisplayName','high moved branch');
    yline(ax,1,'k--','DisplayName','ODE boundary'); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'\beta'); ylabel(ax,'max Re(\lambda[D\Phi])');
    title(ax,sprintf('%s: branches are never connected',pathway(p)));
    legend(ax,'Location','best','FontSize',8);
end
title(layout,'Branch-resolved stability under real pathway tuning');
exportgraphics(fig,fullfile(figureRoot,'06_branch_resolved_stability.pdf'), ...
    'ContentType','vector'); close(fig);
end

function local_plot_sensitivity(data,figureRoot)
fig=figure('Visible','off','Color','w','Position',[100 100 1380 900]);
layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
pathways=["L6","Inhibition"];
for p=1:2
    rows=sortrows(data(data.pathway==pathways(p),:),'beta');
    ax=nexttile(layout,p); hold(ax,'on');
    plot(ax,rows.beta,rows.lambdaPrimeDirect,'-o','LineWidth',1.5, ...
        'DisplayName','direct');
    plot(ax,rows.beta,rows.lambdaPrimeMovement,'-s','LineWidth',1.5, ...
        'DisplayName','movement');
    plot(ax,rows.beta,rows.lambdaPrimeTotal,'-^','LineWidth',1.5, ...
        'DisplayName','sum');
    plot(ax,rows.beta,rows.lambdaPrimeTrackedFiniteDifference,'k--', ...
        'LineWidth',1.3,'DisplayName','tracked finite difference');
    yline(ax,0,'k:'); grid(ax,'on'); xlabel(ax,'\beta');
    ylabel(ax,'d Re(\lambda)/d\beta'); title(ax,pathways(p));
    legend(ax,'Location','best','FontSize',8);
    ax=nexttile(layout,p+2); hold(ax,'on');
    semilogy(ax,rows.beta,rows.directionalResolventGain,'-o','LineWidth',1.5, ...
        'DisplayName','||v||/||b||');
    semilogy(ax,rows.beta,rows.fullResolventTwoNorm,'-s','LineWidth',1.5, ...
        'DisplayName','||(I-J)^{-1}||_2');
    grid(ax,'on'); xlabel(ax,'\beta'); ylabel(ax,'resolvent amplification');
    title(ax,sprintf('%s equilibrium sensitivity',pathways(p)));
    legend(ax,'Location','best','FontSize',8);
end
title(layout,'Direct and operating-point contributions to the tracked eigenvalue');
exportgraphics(fig,fullfile(figureRoot,'07_eigenvalue_sensitivity_decomposition.pdf'), ...
    'ContentType','vector'); close(fig);
end

function local_plot_interpolation(interpolation,shapley,meanPattern,figureRoot)
cases=unique(interpolation(:,{'taskId','pathway','beta'}),'rows','stable');
fig=figure('Visible','off','Color','w','Position',[100 100 1500 1050]);
layout=tiledlayout(fig,3,height(cases),'TileSpacing','compact','Padding','compact');
for caseIndex=1:height(cases)
    mask=interpolation.taskId==cases.taskId(caseIndex);
    rows=sortrows(interpolation(mask,:),'alpha');
    ax=nexttile(layout,caseIndex); plot(ax,rows.alpha,rows.maxRealLambda,'-o', ...
        'LineWidth',1.7); yline(ax,1,'k--'); grid(ax,'on');
    xlabel(ax,'state interpolation \alpha'); ylabel(ax,'max Re(\lambda)');
    title(ax,sprintf('%s \beta=%+.3f',cases.pathway(caseIndex),cases.beta(caseIndex)));
    srows=shapley(shapley.taskId==cases.taskId(caseIndex),:);
    ax=nexttile(layout,height(cases)+caseIndex);
    bar(ax,srows.shapleyContribution); set(ax,'XTick',1:3,'XTickLabel',srows.population);
    yline(ax,0,'k:'); grid(ax,'on'); ylabel(ax,'Shapley contribution to max Re(\lambda)');
    mrows=meanPattern(meanPattern.taskId==cases.taskId(caseIndex),:);
    ax=nexttile(layout,2*height(cases)+caseIndex);
    bar(ax,mrows.maxRealLambda); set(ax,'XTick',1:height(mrows), ...
        'XTickLabel',mrows.stateComponent,'XTickLabelRotation',25);
    yline(ax,1,'k--'); grid(ax,'on'); ylabel(ax,'max Re(\lambda)');
end
title(layout,'Causal state movement: interpolation, population effects, and spatial pattern');
exportgraphics(fig,fullfile(figureRoot,'08_state_movement_causal_diagnostics.pdf'), ...
    'ContentType','vector'); close(fig);
end

function local_plot_threshold(data,figureRoot)
fig=figure('Visible','off','Color','w','Position',[100 100 1450 850]);
layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
pathways=["L6","Inhibition"];
for p=1:2
    subset=data(data.pathway==pathways(p) & data.protocol=="baseline reset",:);
    ax=nexttile(layout,p); hold(ax,'on');
    scales=unique(subset.edgeScale(isfinite(subset.edgeScale)));
    for index=1:numel(scales)
        rows=sortrows(subset(abs(subset.edgeScale-scales(index))<1e-12,:),'beta');
        plot(ax,rows.beta,rows.maximumRateHz,'LineWidth',1.5, ...
            'DisplayName',sprintf('edge x %.1f',scales(index)));
    end
    raw=sortrows(subset(subset.mapMode=="raw",:),'beta');
    plot(ax,raw.beta,raw.maximumRateHz,'k--','LineWidth',1.7,'DisplayName','raw map');
    yline(ax,180,'k:'); grid(ax,'on'); xlabel(ax,'\beta');
    ylabel(ax,'maximum fixed-point rate (Hz)'); title(ax,pathways(p));
    legend(ax,'Location','best','FontSize',8);

    ax=nexttile(layout,p+2); hold(ax,'on');
    for index=1:numel(scales)
        rows=sortrows(subset(abs(subset.edgeScale-scales(index))<1e-12,:),'beta');
        plot(ax,rows.beta,rows.meanExtensionWeight,'LineWidth',1.5, ...
            'DisplayName',sprintf('edge x %.1f',scales(index)));
    end
    grid(ax,'on'); xlabel(ax,'\beta'); ylabel(ax,'mean extension weight');
    title(ax,sprintf('%s extension activation',pathways(p)));
end
title(layout,'Extension-threshold control of the high-rate branch jump');
exportgraphics(fig,fullfile(figureRoot,'09_extension_threshold_control.pdf'), ...
    'ContentType','vector'); close(fig);
end

function local_plot_ode(data,figureRoot)
fig=figure('Visible','off','Color','w','Position',[100 100 900 600]);
ax=axes(fig); hold(ax,'on'); names=unique(data.direction,'stable');
for index=1:numel(names)
    rows=data(data.direction==names(index),:);
    semilogy(ax,rows.timeMs,rows.normalizedLeadingModeAmplitude, ...
        'LineWidth',1.6,'DisplayName',names(index));
end
rows=data(data.direction==names(1),:);
semilogy(ax,rows.timeMs,rows.linearPrediction,'k--','LineWidth',1.5, ...
    'DisplayName','linear prediction');
grid(ax,'on'); xlabel(ax,'ODE time (ms)');
ylabel(ax,'leading left-mode amplitude / initial perturbation');
title(ax,'Nonlinear ODE test of the symmetry-hidden L6 instability');
legend(ax,'Location','northwest','FontSize',8);
exportgraphics(fig,fullfile(figureRoot,'10_hidden_mode_continuous_ode.pdf'), ...
    'ContentType','vector'); close(fig);
end
