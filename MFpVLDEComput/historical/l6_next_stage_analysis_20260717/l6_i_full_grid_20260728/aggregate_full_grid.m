function aggregate_full_grid()
% Assemble all independent grid results and plot canonical HCnorm heatmaps.

outputRoot=getenv('FULL_GRID_OUTPUT_ROOT');
combinedFile=fullfile(outputRoot,'all_beta6_betaI_grid_points.tsv');
if ~isfile(combinedFile)
    error('FullGrid:MissingCombined','Missing combined grid table.');
end
grid=readtable(combinedFile,'FileType','text','Delimiter','\t', ...
    'TextType','string');
gridStep=0.002;
gridMaximum=0.4;
gridCount=round(gridMaximum/gridStep)+1;
totalCount=gridCount^2;
if height(grid)~=totalCount || numel(unique(grid.taskId))~=totalCount
    error('FullGrid:Incomplete','Expected %d unique rows, found %d.', ...
        totalCount,height(grid));
end
grid=sortrows(grid,'taskId');
expectedTask=(1:totalCount)';
if any(grid.taskId~=expectedTask)
    error('FullGrid:Order','Grid task ids are incomplete or out of order.');
end

beta6Values=(0:gridCount-1)*gridStep;
betaIValues=(0:gridCount-1)*gridStep;
hcMatrix=nan(gridCount,gridCount);
convergedMatrix=false(gridCount,gridCount);
residualMatrix=nan(gridCount,gridCount);
iterationMatrix=nan(gridCount,gridCount);
for row=1:height(grid)
    i6=grid.beta6Index(row);
    iI=grid.betaIIndex(row);
    hcMatrix(iI,i6)=grid.HCnorm(row);
    convergedMatrix(iI,i6)=logical(grid.converged(row));
    residualMatrix(iI,i6)=grid.fixedPointResidual(row);
    iterationMatrix(iI,i6)=grid.iterations(row);
end
validHC=hcMatrix;
validHC(~convergedMatrix | ~isfinite(validHC) | validHC<=0)=NaN;

metrics=table( ...
    ["grid_condition_count";"converged_condition_count"; ...
    "nonconverged_condition_count";"minimum_converged_HCnorm"; ...
    "maximum_converged_HCnorm";"maximum_converged_residual"; ...
    "maximum_iterations"], ...
    [totalCount;nnz(convergedMatrix);nnz(~convergedMatrix); ...
    min(validHC,[],'all','omitnan');max(validHC,[],'all','omitnan'); ...
    max(residualMatrix(convergedMatrix),[],'omitnan'); ...
    max(iterationMatrix,[],'all','omitnan')], ...
    'VariableNames',{'metric','value'});
writetable(metrics,fullfile(outputRoot,'full_grid_metrics.tsv'), ...
    'FileType','text','Delimiter','\t');

figureRoot=fullfile(outputRoot,'figures');
if ~exist(figureRoot,'dir'); mkdir(figureRoot); end
theorySlope=str2double(getenv('FULL_GRID_THEORY_SLOPE'));
if ~isfinite(theorySlope); theorySlope=1.30243705978516; end

fig=figure('Visible','off','Color','w','Position',[80 80 1280 900]);
heatmapImage=imagesc(beta6Values,betaIValues,validHC);
heatmapImage.AlphaData=isfinite(validHC);
set(gca,'YDir','normal','Color',[0.80 0.80 0.80]);
axis xy tight; hold on; set(gca,'XGrid','off','YGrid','off');
lineX=linspace(0,min(gridMaximum,gridMaximum/theorySlope),300);
plot(lineX,theorySlope*lineX,'k--','LineWidth',2.2, ...
    'DisplayName',sprintf('theoretical local stability slope %.6f',theorySlope));
xlabel('\beta_6 (fractional whole-L6 increase)', ...
    'FontSize',18,'FontWeight','bold');
ylabel('\beta_I (fractional whole-inhibition increase)', ...
    'FontSize',18,'FontWeight','bold');
set(gca,'FontSize',16,'LineWidth',1.1);
colormap(fig,jet);
ax=gca; ax.ColorScale='log';
finiteValues=validHC(isfinite(validHC));
finiteValues=sort(finiteValues);
lowerIndex=max(1,round(0.001*(numel(finiteValues)-1))+1);
upperIndex=min(numel(finiteValues),round(0.98*(numel(finiteValues)-1))+1);
clim([finiteValues(lowerIndex) finiteValues(upperIndex)]);
cb=colorbar;
cb.Label.String='HC norm';
cb.Label.FontSize=16;
cb.Label.FontWeight='bold';
cb.Title.String='Hz';
cb.Title.FontSize=15;
cb.Title.FontWeight='bold';
legend('Location','northwest','Interpreter','none');
local_save(fig,figureRoot,'07_FullGrid_Beta6_BetaI_Canonical_HCnorm_Heatmap');

save(fullfile(outputRoot,'full_grid_results.mat'),'grid','beta6Values', ...
    'betaIValues','hcMatrix','validHC','convergedMatrix','residualMatrix', ...
    'iterationMatrix','metrics','theorySlope','-v7.3');
fprintf('Aggregated %d grid points: %d converged, %d nonconverged.\n', ...
    totalCount,nnz(convergedMatrix),nnz(~convergedMatrix));
end

function local_save(fig,root,stem)
savefig(fig,fullfile(root,[stem '.fig']));
exportgraphics(fig,fullfile(root,[stem '.pdf']),'ContentType','vector');
close(fig);
end
