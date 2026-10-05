function result = analyze_l6ns_inhibition_active_set_event(root)
% Document the clamp transition replacing the former smooth-fold claim.

if nargin<1 || isempty(root); root=pwd; end
loaded=load(fullfile(root,'branch_atlas_20260803','smoke', ...
    'inhibition_fold_ultrafine.mat'),'result');
branch=loaded.result;
loaded=load(fullfile(root,'results_global_bifurcation_20260722', ...
    'global_bifurcation_setup.mat'));
names=fieldnames(loaded); setup=loaded.(names{1});
fraction=branch.Summary.clampActiveFraction;
transition=find(abs(diff(fraction))>0,1,'first');
if isempty(transition); error('No L6 clamp active-set transition was found.'); end
[~,rawBefore]=l6ns_l6_dynamic(branch.States{transition},setup.Context,[1 1]);
[~,rawAfter]=l6ns_l6_dynamic(branch.States{transition+1},setup.Context,[1 1]);
activeBefore=rawBefore<=1 | rawBefore>=40;
activeAfter=rawAfter<=1 | rawAfter>=40;
enter=find(~activeBefore & activeAfter);
leave=find(activeBefore & ~activeAfter);
[leaveRow,leaveColumn]=ind2sub(setup.Context.MapSize,leave);
siteTable=table(leave(:),leaveRow(:),leaveColumn(:),rawBefore(leave),rawAfter(leave), ...
    'VariableNames',{'linearIndex','row','column','rawL6Before','rawL6After'});
summary=table(branch.Summary.freezeWeight(transition), ...
    branch.Summary.freezeWeight(transition+1), ...
    branch.Summary.maxRealLambda(transition), ...
    branch.Summary.maxRealLambda(transition+1), ...
    branch.Summary.unstableDimension(transition), ...
    branch.Summary.unstableDimension(transition+1), ...
    fraction(transition),fraction(transition+1),numel(enter),numel(leave), ...
    string('piecewise_smooth_lower_L6_clamp_border_event'), ...
    'VariableNames',{'weightBefore','weightAfter','maxRealBefore','maxRealAfter', ...
    'unstableDimensionBefore','unstableDimensionAfter','activeFractionBefore', ...
    'activeFractionAfter','sitesEnteringClamp','sitesLeavingClamp','classification'});
result=struct('Summary',summary,'SitesLeavingClamp',siteTable, ...
    'BeforeState',branch.States{transition},'AfterState',branch.States{transition+1});
outputDir=fullfile(root,'branch_atlas_20260803');
save(fullfile(outputDir,'inhibition_active_set_event.mat'),'result','-v7.3');
writetable(summary,fullfile(outputDir,'inhibition_active_set_event.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(siteTable,fullfile(outputDir,'inhibition_active_set_event_sites.tsv'), ...
    'FileType','text','Delimiter','\t');

figureDir=fullfile(outputDir,'figures');
if ~exist(figureDir,'dir'); mkdir(figureDir); end
figureHandle=figure('Color','w','Position',[80 80 1450 850]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; imagesc(reshape(rawBefore,setup.Context.MapSize)); axis image;
clim([0.9995 1.0005]); colorbar; title('Raw L6 immediately before release');
xlabel('column'); ylabel('row');
nexttile; imagesc(reshape(rawAfter,setup.Context.MapSize)); axis image;
clim([0.9995 1.0005]); colorbar; title('Raw L6 immediately after release');
xlabel('column'); ylabel('row');
nexttile; imagesc(reshape(activeBefore & ~activeAfter,setup.Context.MapSize));
axis image; colorbar; title('32 sites leaving the lower clamp');
xlabel('column'); ylabel('row');
nexttile;
window=max(1,transition-5):min(height(branch.Summary),transition+6);
yyaxis left;
plot(branch.Summary.freezeWeight(window),branch.Summary.maxRealLambda(window), ...
    'o-','LineWidth',1.4); hold on; yline(1,'r--');
ylabel('max Re lambda(J)');
yyaxis right;
stairs(branch.Summary.freezeWeight(window),branch.Summary.clampActiveFraction(window), ...
    's-','LineWidth',1.4);
ylabel('L6 clamp-active fraction'); xlabel('inhibition freeze weight w_I');
title('Jacobian jump at the active-set boundary'); grid on; box on;
sgtitle('Inhibition near-silent branch: piecewise-smooth L6 clamp event');
stem=fullfile(figureDir,'inhibition_nearsilent_active_set_border_event');
exportgraphics(figureHandle,[stem '.pdf'],'ContentType','vector');
savefig(figureHandle,[stem '.fig']);
end
