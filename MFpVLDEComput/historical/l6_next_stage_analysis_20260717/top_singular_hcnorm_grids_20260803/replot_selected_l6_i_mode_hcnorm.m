function replot_selected_l6_i_mode_hcnorm()
% Replot the selected-condition mode-HC map from the saved grid only.

outputRoot = fullfile( ...
    '/Users/xialeihuang/Desktop/Neuroscience_Project/Spectral_Analysis', ...
    'matlab-inserting_into_CG_model/Complete_Code_for_Paper3', ...
    'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors', ...
    'L6 and Inhibition/L6_and_Inhibition_offset');
stem = '6.1.4a_Selected_Conditions_on_Mode_HCnorm_Heatmap';
points = readtable(fullfile(outputRoot,'all_top_singular_grid_points.tsv'), ...
    'FileType','text','Delimiter','\t');
gridValues = 0:0.008:0.4;
modeHC = nan(numel(gridValues));
for index = 1:height(points)
    if points.fixedPointAccepted(index) && isfinite(points.modeHCnorm(index))
        modeHC(points.betaYIndex(index),points.betaXIndex(index)) = ...
            points.modeHCnorm(index);
    end
end

selectedBeta6 = [0.024 0.336 0.080 0.272 0.352 0.376];
selectedBetaI = [0.376 0.376 0.024 0.016 0.056 0.008];
labels = ["A" "B" "C" "D" "E" "F"];

fig = figure('Visible','off','Color','w','Position',[80 80 1280 1020]);
imageHandle = imagesc(gridValues,gridValues,modeHC);
imageHandle.AlphaData = isfinite(modeHC);
ax = gca;
set(ax,'YDir','normal','Color',[0.80 0.80 0.80]);
axis(ax,'xy');
axis(ax,'tight');
xlabel(ax,'\beta_6 (L6 increase)','Interpreter','tex', ...
    'FontSize',22,'FontWeight','bold');
ylabel(ax,'\beta_I (inhibition increase)','Interpreter','tex', ...
    'FontSize',22,'FontWeight','bold');
set(ax,'FontSize',19,'FontWeight','bold','LineWidth',1.25);
colormap(fig,jet(256));
cb = colorbar(ax);
cb.FontSize = 17;
cb.Label.String = 'Mode HC norm';
cb.Label.FontSize = 19;
cb.Label.FontWeight = 'bold';
hold(ax,'on');

bestBetaI = nan(size(gridValues));
for xIndex = 1:numel(gridValues)
    columnHC = modeHC(:,xIndex);
    finiteRows = find(isfinite(columnHC));
    if isempty(finiteRows); continue; end
    [~,localIndex] = min(columnHC(finiteRows));
    bestBetaI(xIndex) = gridValues(finiteRows(localIndex));
end
fitMask = isfinite(bestBetaI);
regressionSlope = sum(gridValues(fitMask).*bestBetaI(fitMask)) / ...
    sum(gridValues(fitMask).^2);
lineX = linspace(min(gridValues),max(gridValues),500);
lineY = regressionSlope*lineX;
inside = lineY >= min(gridValues) & lineY <= max(gridValues);
hRegression = plot(ax,lineX(inside),lineY(inside),'-', ...
    'Color',[0.35 0.35 0.35],'LineWidth',2.2);

for index = 1:numel(selectedBeta6)
    plot(ax,selectedBeta6(index),selectedBetaI(index),'s', ...
        'MarkerSize',15,'LineWidth',5,'Color','w');
    plot(ax,selectedBeta6(index),selectedBetaI(index),'s', ...
        'MarkerSize',15,'LineWidth',2.5,'Color','k');
    text(ax,selectedBeta6(index)+0.008,selectedBetaI(index),labels(index), ...
        'Color','k','BackgroundColor','w','Margin',1, ...
        'FontSize',15,'FontWeight','bold','VerticalAlignment','middle');
end
legend(ax,hRegression,sprintf('minimum-mode-HC regression %.6f', ...
    regressionSlope),'Location','northoutside','FontSize',12, ...
    'FontWeight','bold','Box','off');
savefig(fig,fullfile(outputRoot,[stem '.fig']));
exportgraphics(fig,fullfile(outputRoot,[stem '.pdf']), ...
    'ContentType','image','Resolution',300);
close(fig);

summary = table(regressionSlope, ...
    'VariableNames',{'minimumModeHCRegressionSlope'});
writetable(summary,fullfile(outputRoot, ...
    '6.1.4a_Mode_HCnorm_Regression.tsv'),'FileType','text', ...
    'Delimiter','\t');
disp(summary);
end
