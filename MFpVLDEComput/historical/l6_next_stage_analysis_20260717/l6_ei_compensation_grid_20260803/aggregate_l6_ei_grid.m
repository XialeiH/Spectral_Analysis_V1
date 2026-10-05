function aggregate_l6_ei_grid()
% Assemble the 41-by-71 betaEI/beta6 scan and generate Figure 6.1.2.

outputRoot=getenv('L6EI_OUTPUT_ROOT');
setupFile=getenv('L6EI_SETUP_FILE');
combinedFile=fullfile(outputRoot,'all_beta6_betaEI_grid_points.tsv');
if ~isfile(combinedFile)
    error('L6EI:MissingCombined','Missing combined grid table.');
end
grid=readtable(combinedFile,'FileType','text','Delimiter','\t', ...
    'TextType','string');
gridStep=local_env_number('L6EI_GRID_STEP',0.005);
beta6Maximum=local_env_number('L6EI_BETA6_MAXIMUM',0.35);
betaEIMaximum=local_env_number('L6EI_BETAEI_MAXIMUM',0.2);
beta6Count=round(beta6Maximum/gridStep)+1;
betaEICount=round(betaEIMaximum/gridStep)+1;
totalCount=beta6Count*betaEICount;
if height(grid)~=totalCount || numel(unique(grid.taskId))~=totalCount
    error('L6EI:Incomplete','Expected %d unique rows, found %d.', ...
        totalCount,height(grid));
end
grid=sortrows(grid,'taskId');
if any(grid.taskId~=(1:totalCount)')
    error('L6EI:Order','Grid task ids are incomplete or out of order.');
end

beta6Values=(0:beta6Count-1)*gridStep;
betaEIValues=(0:betaEICount-1)*gridStep;
hcMatrix=nan(betaEICount,beta6Count);
acceptedMatrix=false(betaEICount,beta6Count);
strictMatrix=false(betaEICount,beta6Count);
residualMatrix=nan(betaEICount,beta6Count);
iterationMatrix=nan(betaEICount,beta6Count);
for row=1:height(grid)
    i6=grid.beta6Index(row);
    iEI=grid.betaEIIndex(row);
    hcMatrix(iEI,i6)=grid.HCnorm(row);
    acceptedMatrix(iEI,i6)=logical(grid.accepted(row));
    strictMatrix(iEI,i6)=logical(grid.strictConverged(row));
    residualMatrix(iEI,i6)=grid.fixedPointResidual(row);
    iterationMatrix(iEI,i6)=grid.iterations(row);
end
validHC=hcMatrix;
validHC(~acceptedMatrix | ~isfinite(validHC) | validHC<0)=NaN;
positiveHC=validHC(isfinite(validHC) & validHC>0);
if isempty(positiveHC)
    error('L6EI:NoPositiveHC','No positive accepted HCnorm values exist.');
end
plotHC=validHC;
plotHC(isfinite(plotHC) & plotHC<=0)=realmin('double');

bestBetaEI=nan(beta6Count,1);
bestHC=nan(beta6Count,1);
for i6=1:beta6Count
    column=validHC(:,i6);
    [bestHC(i6),iEI]=min(column,[],'omitnan');
    if isfinite(bestHC(i6)); bestBetaEI(i6)=betaEIValues(iEI); end
end
bestPath=table(beta6Values(:),bestBetaEI,bestHC, ...
    'VariableNames',{'beta6','betaEI','minimumHCnorm'});
writetable(bestPath,fullfile(outputRoot,'best_compensation_path.tsv'), ...
    'FileType','text','Delimiter','\t');

audit=l6ei_first_order_slope(setupFile,outputRoot);
theorySlope=audit.slopeL6ItoE(1);
[minimumHC,linearMinimumIndex]=min(validHC(:),[],'omitnan');
[minimumEIIndex,minimum6Index]=ind2sub(size(validHC),linearMinimumIndex);
metrics=table( ...
    ["grid_condition_count";"strict_converged_condition_count"; ...
    "accepted_condition_count";"gray_condition_count"; ...
    "accepted_HCnorm_le_3_count";"minimum_accepted_HCnorm"; ...
    "minimum_HC_beta6";"minimum_HC_betaEI"; ...
    "maximum_accepted_HCnorm";"maximum_accepted_residual"; ...
    "maximum_iterations";"theoretical_FPP_slope_betaEI_per_beta6"], ...
    [totalCount;nnz(strictMatrix);nnz(acceptedMatrix); ...
    nnz(~acceptedMatrix);nnz(validHC<=3);minimumHC; ...
    beta6Values(minimum6Index);betaEIValues(minimumEIIndex); ...
    max(validHC,[],'all','omitnan'); ...
    max(residualMatrix(acceptedMatrix),[],'omitnan'); ...
    max(iterationMatrix,[],'all','omitnan');theorySlope], ...
    'VariableNames',{'metric','value'});
writetable(metrics,fullfile(outputRoot,'full_grid_metrics.tsv'), ...
    'FileType','text','Delimiter','\t');

figureRoot=fullfile(outputRoot,'figures');
if ~exist(figureRoot,'dir'); mkdir(figureRoot); end
fig=figure('Visible','off','Color','w','Position',[80 80 1280 900]);
heatmapImage=imagesc(beta6Values,betaEIValues,plotHC);
heatmapImage.AlphaData=isfinite(plotHC);
set(gca,'YDir','normal','Color',[0.80 0.80 0.80]);
axis xy tight
hold on
if min(validHC,[],'all','omitnan')<=3 && ...
        max(validHC,[],'all','omitnan')>=3
    contour(beta6Values,betaEIValues,validHC,[3 3], ...
        'k-','LineWidth',1.25,'DisplayName','HC norm = 3 Hz');
end
lineX=linspace(0,min(beta6Maximum,betaEIMaximum/theorySlope),300);
plot(lineX,theorySlope*lineX,'k--','LineWidth',2.2, ...
    'DisplayName',sprintf('first-order FPP slope %.3f',theorySlope));
xlabel('\beta_6 (fractional whole-L6 increase)', ...
    'FontSize',22,'FontWeight','bold');
ylabel('\beta_{EI} (fractional I-to-E inhibition increase)', ...
    'FontSize',22,'FontWeight','bold');
set(gca,'FontSize',19,'FontWeight','bold','LineWidth',1.25, ...
    'XGrid','off','YGrid','off');
colormap(fig,jet(256));
ax=gca;
ax.ColorScale='log';
referenceClim=[0.391523346338067 36.6263881897836];
clim(referenceClim);
cb=colorbar;
cb.FontSize=17;
cb.Label.String='HC norm';
cb.Label.FontSize=19;
cb.Label.FontWeight='bold';
cb.Title.String='Hz';
cb.Title.FontSize=17;
cb.Title.FontWeight='bold';
legend('Location','northoutside','Orientation','horizontal', ...
    'FontSize',11.5,'Box','off');
stem='6.1.2_FullGrid_Beta6_BetaEI_Canonical_HCnorm_Heatmap_with_HCnorm3_Boundary';
savefig(fig,fullfile(figureRoot,[stem '.fig']));
exportgraphics(fig,fullfile(figureRoot,[stem '.pdf']), ...
    'ContentType','image','Resolution',300);
close(fig);

save(fullfile(outputRoot,'full_grid_results.mat'),'grid','beta6Values', ...
    'betaEIValues','hcMatrix','validHC','plotHC','acceptedMatrix', ...
    'strictMatrix','residualMatrix','iterationMatrix','bestPath', ...
    'metrics','audit','theorySlope','referenceClim','-v7.3');
fprintf(['Aggregated %d L6/I-to-E points: %d accepted at residual <=1e-5; ' ...
    '%d gray.\n'],totalCount,nnz(acceptedMatrix),nnz(~acceptedMatrix));
end

function value=local_env_number(name,defaultValue)
value=str2double(getenv(name));
if ~isfinite(value); value=defaultValue; end
end
