function aggregate_l6_i_offset()
% Aggregate verified offsets, regress beta_I on beta_6, and plot results.

outputRoot = getenv('OFFSET_OUTPUT_ROOT');
setupFile = getenv('OFFSET_SETUP_FILE');
if isempty(outputRoot) || ~exist(outputRoot,'dir') || ~isfile(setupFile)
    error('Offset:Environment','OFFSET_OUTPUT_ROOT and OFFSET_SETUP_FILE are required.');
end
files = dir(fullfile(outputRoot,'*','case_summary.tsv'));
if isempty(files); error('Offset:NoResults','No case summaries found.'); end
parts = cell(numel(files),1);
for index=1:numel(files)
    parts{index}=readtable(fullfile(files(index).folder,files(index).name), ...
        'FileType','text','Delimiter','\t','TextType','string');
end
summary = sortrows(vertcat(parts{:}),'beta6');
writetable(summary,fullfile(outputRoot,'offset_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

valid = summary.converged & summary.HCnorm<0.3;
if nnz(valid)<3; error('Offset:TooFewPoints','Need at least three HC-valid points.'); end
x = summary.beta6(valid);
y = summary.betaI(valid);
ordinary = polyfit(x,y,1);
yOrdinary = polyval(ordinary,x);
r2 = 1-sum((y-yOrdinary).^2)/sum((y-mean(y)).^2);
positive = x>0;
slopeOrigin = dot(x(positive),y(positive))/dot(x(positive),x(positive));
residualOrigin = y(positive)-slopeOrigin*x(positive);
r2Origin = 1-sum(residualOrigin.^2)/sum(y(positive).^2);
slopeOriginSE = sqrt(sum(residualOrigin.^2)/max(nnz(positive)-1,1) / ...
    dot(x(positive),x(positive)));

metrics = table(["through_origin_slope";"through_origin_slope_se"; ...
    "through_origin_r2";"ordinary_slope";"ordinary_intercept"; ...
    "ordinary_r2";"valid_point_count";"maximum_valid_beta6"], ...
    [slopeOrigin;slopeOriginSE;r2Origin;ordinary(1);ordinary(2);r2; ...
    nnz(valid);max(x)],'VariableNames',{'metric','value'});
writetable(metrics,fullfile(outputRoot,'offset_regression.tsv'), ...
    'FileType','text','Delimiter','\t');

loaded=load(setupFile,'setup');
contrast=unique(loaded.setup.Context.ContrastUse(:));
angle=unique(loaded.setup.Context.OrientationUse(:));
condition=sprintf('angle %.2f deg, contrast %.4g',angle(1),contrast(1));
figureRoot=fullfile(outputRoot,'figures');
if ~exist(figureRoot,'dir'); mkdir(figureRoot); end

fig=figure('Visible','off','Color','w','Position',[100 100 1050 700]);
hold on; grid on;
plot(summary.beta6(~valid),summary.betaI(~valid),'x','Color',[0.45 0.45 0.45], ...
    'LineWidth',1.8,'MarkerSize',8,'DisplayName','optimized but HCnorm >= 0.3');
scatter(x,y,55,[0.08 0.42 0.72],'filled','DisplayName','HCnorm < 0.3');
xLine=linspace(0,max(summary.beta6),250);
plot(xLine,slopeOrigin*xLine,'-','Color',[0.85 0.33 0.10],'LineWidth',2.2, ...
    'DisplayName',sprintf('through origin: betaI = %.4f beta6',slopeOrigin));
plot(xLine,polyval(ordinary,xLine),'--','Color',[0.20 0.60 0.25], ...
    'LineWidth',1.8,'DisplayName',sprintf('ordinary: betaI = %.4f beta6 %+.2g', ...
    ordinary(1),ordinary(2)));
xlabel('\beta_6 (fractional whole-L6 increase)');
ylabel('\beta_I (fractional whole-inhibition increase)');
title({sprintf('Real L6-inhibition offset with moved fixed point: %s',condition), ...
    sprintf(['Unmoved criterion: HCnorm < 0.3 | through-origin slope ' ...
    '%.5f | R2 %.5f'],slopeOrigin,r2Origin)}, ...
    'FontWeight','bold','Interpreter','none');
legend('Location','northwest','Interpreter','none');
save_figure(fig,figureRoot,'01_L6_Inhibition_Offset_Regression');

fig=figure('Visible','off','Color','w','Position',[100 100 1050 700]);
plot(summary.beta6,summary.HCnorm,'-o','Color',[0.08 0.42 0.72], ...
    'LineWidth',2,'MarkerSize',5); hold on; grid on;
yline(0.3,'r--','LineWidth',1.8,'Label','unmoved threshold');
xlabel('\beta_6 (fractional whole-L6 increase)'); ylabel('minimum verified HCnorm');
title({sprintf('Best fixed-point cancellation after optimizing betaI: %s',condition), ...
    sprintf('Largest sampled beta6 satisfying HCnorm < 0.3: %.4f',max(x))}, ...
    'FontWeight','bold','Interpreter','none');
save_figure(fig,figureRoot,'02_Minimum_HCnorm_vs_Beta6');

save(fullfile(outputRoot,'offset_results.mat'),'summary','valid','metrics', ...
    'ordinary','r2','slopeOrigin','slopeOriginSE','r2Origin','condition','-v7.3');
fprintf(['Aggregated %d offsets: %d satisfy HC<0.3; origin slope=%.8f ' ...
    '(SE %.3g), ordinary slope=%.8f intercept=%+.3g.\n'],height(summary), ...
    nnz(valid),slopeOrigin,slopeOriginSE,ordinary(1),ordinary(2));
end

function save_figure(fig,root,stem)
savefig(fig,fullfile(root,[stem '.fig']));
exportgraphics(fig,fullfile(root,[stem '.pdf']),'ContentType','vector');
close(fig);
end
