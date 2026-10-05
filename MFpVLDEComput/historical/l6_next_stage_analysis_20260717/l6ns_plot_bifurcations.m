function l6ns_plot_bifurcations(result,outputDir)
% Plot broad stability and magnified local branches for L6 and inhibition.

local_plot_one(result.L6,outputDir,'l6_bifurcation_diagram', ...
    'L6 dynamic feedback','w_6');
local_plot_one(result.Inhibition,outputDir,'inhibition_bifurcation_diagram', ...
    'inhibitory-output feedback','w_I');
local_plot_combined(result.L6,result.Inhibition,outputDir);
end

function local_plot_one(result,outputDir,stem,pathwayTitle,weightLabel)
fig = figure('Visible','off','Color','w','Position',[80 60 1150 900]);
layout = tiledlayout(2,1,'TileSpacing','loose','Padding','loose');

ax1 = nexttile(layout);
curve = result.StabilityCurve;
plot(ax1,curve.freezeWeight,curve.maxRealLambda,'k-','LineWidth',1.8); hold(ax1,'on');
yline(ax1,1,'--','Color',[0.75 0.1 0.1],'LineWidth',1.2);
xline(ax1,result.CriticalFreezeWeight,':','Color',[0.2 0.4 0.8],'LineWidth',1.4);
grid(ax1,'on');
xlabel(ax1,sprintf('%s (freeze weight)',weightLabel));
ylabel(ax1,'max Re \lambda(D\Phi)');
title(ax1,sprintf('%s: persistent fixed-point branch, angle 0 deg, contrast 100', ...
    pathwayTitle),'FontWeight','bold');
legend(ax1,{'persistent branch','stability boundary','critical weight'}, ...
    'Location','best');

ax2 = nexttile(layout);
local_plot_branch_axis(ax2,result,weightLabel,true);
title(ax2,sprintf('Magnified local branches: %s', ...
    strrep(result.Classification,'_',' ')),'FontWeight','bold');
set([ax1 ax2],'FontSize',11);
l6ns_save_figure(fig,outputDir,stem); close(fig);
end

function local_plot_combined(l6,inhibition,outputDir)
fig = figure('Visible','off','Color','w','Position',[50 50 1400 940]);
layout = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
local_plot_stability_axis(nexttile(layout),l6,'L6','w_6');
local_plot_branch_axis(nexttile(layout),l6,'w_6',false);
local_plot_stability_axis(nexttile(layout),inhibition,'inhibition','w_I');
local_plot_branch_axis(nexttile(layout),inhibition,'w_I',false);
sgtitle(layout,'Fixed-point-preserving L6 and inhibition bifurcation comparison', ...
    'FontWeight','bold','FontSize',15);
l6ns_save_figure(fig,outputDir,'l6_and_inhibition_bifurcation_comparison'); close(fig);
end

function local_plot_stability_axis(ax,result,label,weightLabel)
plot(ax,result.StabilityCurve.freezeWeight,result.StabilityCurve.maxRealLambda, ...
    'k-','LineWidth',1.7); hold(ax,'on');
yline(ax,1,'r--','LineWidth',1.1);
xline(ax,result.CriticalFreezeWeight,'b:','LineWidth',1.2);
grid(ax,'on'); xlabel(ax,weightLabel); ylabel(ax,'max Re \lambda(D\Phi)');
title(ax,sprintf('%s persistent branch',label),'FontWeight','bold');
set(ax,'FontSize',11);
end

function local_plot_branch_axis(ax,result,weightLabel,showLegend)
hold(ax,'on');
t = result.Continuation;
valid = t.converged & isfinite(t.maxRealLambda);
xBranch = 1e6*(t.freezeWeight(valid)-result.CriticalFreezeWeight);
scale = max(abs(t.amplitude(valid)));
if isempty(scale) || scale==0; scale=1; end
[xMinimum,xMaximum] = local_branch_limits(xBranch);

belowStable = local_side_stability(result,-1);
aboveStable = local_side_stability(result,+1);
local_persistent_segment(ax,[xMinimum 0],belowStable);
local_persistent_segment(ax,[0 xMaximum],aboveStable);

local_secondary(ax,t,valid & t.amplitude<0,result.CriticalFreezeWeight, ...
    scale,[0.10 0.45 0.80],'secondary a<0');
local_secondary(ax,t,valid & t.amplitude>0,result.CriticalFreezeWeight, ...
    scale,[0.85 0.33 0.10],'secondary a>0');
xline(ax,0,':','Color',[0.2 0.4 0.8],'LineWidth',1.2, ...
    'DisplayName','critical weight');
xlim(ax,[xMinimum xMaximum]); ylim(ax,[-1.08 1.08]);
grid(ax,'on');
xlabel(ax,sprintf('10^6 (%s-%s^c)',weightLabel,weightLabel));
ylabel(ax,'signed critical amplitude / max |a|');
if showLegend
    legend(ax,'show','Location','best');
else
    if strcmp(weightLabel,'w_6'); branchTitle='L6'; else; branchTitle='inhibition'; end
    title(ax,sprintf('%s local secondary branches',branchTitle),'FontWeight','bold');
end
set(ax,'FontSize',11);
end

function [xMinimum,xMaximum] = local_branch_limits(x)
xMinimum = min([0;x(:)]);
xMaximum = max([0;x(:)]);
span = max(xMaximum-xMinimum,0.1);
if xMinimum>=-eps
    xMinimum = -0.30*max(xMaximum,span);
elseif xMaximum<=eps
    xMaximum = 0.30*max(abs(xMinimum),span);
end
padding = 0.10*(xMaximum-xMinimum);
xMinimum = xMinimum-padding;
xMaximum = xMaximum+padding;
end

function stable = local_side_stability(result,side)
w = result.StabilityCurve.freezeWeight;
alpha = result.StabilityCurve.maxRealLambda;
if side<0
    candidates = find(w<result.CriticalFreezeWeight);
    [~,index] = max(w(candidates));
else
    candidates = find(w>result.CriticalFreezeWeight);
    [~,index] = min(w(candidates));
end
stable = alpha(candidates(index))<1;
end

function local_persistent_segment(ax,x,stable)
if stable
    style = '-'; label = 'persistent stable';
else
    style = '--'; label = 'persistent unstable';
end
plot(ax,x,[0 0],style,'Color',[0 0 0],'LineWidth',2.2,'DisplayName',label);
end

function local_secondary(ax,t,mask,criticalWeight,scale,color,label)
subset = t(mask,:);
[~,order] = sort(subset.freezeWeight);
subset = subset(order,:);
x = 1e6*(subset.freezeWeight-criticalWeight);
y = subset.amplitude/scale;
if all(subset.maxRealLambda<1); style='-o'; else; style='--o'; end
plot(ax,x,y,style,'Color',color,'MarkerFaceColor',color,'LineWidth',1.7, ...
    'DisplayName',label);
end
