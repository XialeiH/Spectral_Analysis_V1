function aggregate_pair_grid()
% Assemble one 81-by-81 pathway-pair scan and plot canonical HCnorm.

outputRoot=getenv('PAIR_GRID_OUTPUT_ROOT');
pairType=string(getenv('PAIR_GRID_PAIR'));
combinedFile=fullfile(outputRoot,'all_pair_grid_points.tsv');
if ~isfile(combinedFile)
    error('PairGrid:MissingCombined','Missing combined grid table.');
end
gridData=readtable(combinedFile,'FileType','text','Delimiter','\t', ...
    'TextType','string');
gridStep=0.005;
gridMaximum=0.4;
gridCount=round(gridMaximum/gridStep)+1;
totalCount=gridCount^2;
if height(gridData)~=totalCount || ...
        numel(unique(gridData.taskId))~=totalCount
    error('PairGrid:Incomplete','Expected %d unique rows, found %d.', ...
        totalCount,height(gridData));
end
gridData=sortrows(gridData,'taskId');
if any(gridData.taskId~=(1:totalCount)')
    error('PairGrid:Order','Grid task ids are incomplete or out of order.');
end

switch pairType
    case "l6_excitation"
        betaXName='beta6';
        betaYName='betaE';
        xLabel='\beta_6 (fractional whole-L6 increase)';
        yLabel='\beta_E (fractional recurrent-excitation change)';
        stem='08_FullGrid_Beta6_BetaE_Canonical_HCnorm_Heatmap';
        betaXValues=(0:gridCount-1)*gridStep;
        betaYValues=-gridMaximum+(0:gridCount-1)*gridStep;
    case "excitation_inhibition"
        betaXName='betaE';
        betaYName='betaI';
        xLabel='\beta_E (fractional recurrent-excitation increase)';
        yLabel='\beta_I (fractional whole-inhibition increase)';
        stem='09_FullGrid_BetaE_BetaI_Canonical_HCnorm_Heatmap';
        betaXValues=(0:gridCount-1)*gridStep;
        betaYValues=(0:gridCount-1)*gridStep;
    otherwise
        error('PairGrid:PairType','Unsupported pair type: %s.',pairType);
end
hcMatrix=nan(gridCount,gridCount);
convergedMatrix=false(gridCount,gridCount);
residualMatrix=nan(gridCount,gridCount);
iterationMatrix=nan(gridCount,gridCount);
for row=1:height(gridData)
    ix=gridData.betaXIndex(row);
    iy=gridData.betaYIndex(row);
    hcMatrix(iy,ix)=gridData.HCnorm(row);
    convergedMatrix(iy,ix)=logical(gridData.converged(row));
    residualMatrix(iy,ix)=gridData.fixedPointResidual(row);
    iterationMatrix(iy,ix)=gridData.iterations(row);
end
validHC=hcMatrix;
validHC(~convergedMatrix | ~isfinite(validHC) | validHC<0)=NaN;
positiveHC=validHC(isfinite(validHC) & validHC>0);
if isempty(positiveHC)
    error('PairGrid:NoPositiveHC','No positive converged HCnorm values exist.');
end
plotHC=validHC;
zeroFloor=max(min(positiveHC)/2,realmin('double'));
plotHC(plotHC==0)=zeroFloor;

bestBetaY=nan(gridCount,1);
bestHC=nan(gridCount,1);
for ix=1:gridCount
    column=validHC(:,ix);
    [bestHC(ix),iy]=min(column,[],'omitnan');
    if isfinite(bestHC(ix)); bestBetaY(ix)=betaYValues(iy); end
end
bestPath=table(betaXValues(:),bestBetaY,bestHC, ...
    'VariableNames',{betaXName,betaYName,'minimumHCnorm'});
writetable(bestPath,fullfile(outputRoot,'best_compensation_path.tsv'), ...
    'FileType','text','Delimiter','\t');

[minimumHC,linearMinimumIndex]=min(validHC(:),[],'omitnan');
[minimumYIndex,minimumXIndex]=ind2sub(size(validHC),linearMinimumIndex);
metrics=table( ...
    ["grid_condition_count";"converged_condition_count"; ...
    "nonconverged_condition_count";"minimum_converged_HCnorm"; ...
    "minimum_HC_betaX";"minimum_HC_betaY"; ...
    "maximum_converged_HCnorm";"maximum_converged_residual"; ...
    "maximum_iterations"], ...
    [totalCount;nnz(convergedMatrix);nnz(~convergedMatrix);minimumHC; ...
    betaXValues(minimumXIndex);betaYValues(minimumYIndex); ...
    max(validHC,[],'all','omitnan'); ...
    max(residualMatrix(convergedMatrix),[],'omitnan'); ...
    max(iterationMatrix,[],'all','omitnan')], ...
    'VariableNames',{'metric','value'});
writetable(metrics,fullfile(outputRoot,'full_grid_metrics.tsv'), ...
    'FileType','text','Delimiter','\t');

figureRoot=fullfile(outputRoot,'figures');
if ~exist(figureRoot,'dir'); mkdir(figureRoot); end
fig=figure('Visible','off','Color','w','Position',[80 80 1280 1020]);
heatmapImage=imagesc(betaXValues,betaYValues,plotHC);
heatmapImage.AlphaData=isfinite(plotHC);
set(gca,'YDir','normal','Color',[0.80 0.80 0.80]);
axis xy tight
xlabel(xLabel,'FontSize',22,'FontWeight','bold');
ylabel(yLabel,'FontSize',22,'FontWeight','bold');
set(gca,'FontSize',19,'FontWeight','bold','LineWidth',1.25, ...
    'XGrid','off','YGrid','off');
colormap(fig,jet);
ax=gca;
ax.ColorScale='log';
finiteValues=sort(plotHC(isfinite(plotHC)));
lowerIndex=max(1,round(0.001*(numel(finiteValues)-1))+1);
upperIndex=min(numel(finiteValues),round(0.98*(numel(finiteValues)-1))+1);
if finiteValues(lowerIndex)<finiteValues(upperIndex)
    clim([finiteValues(lowerIndex) finiteValues(upperIndex)]);
end
cb=colorbar;
cb.FontSize=17;
cb.Label.String='HC norm';
cb.Label.FontSize=19;
cb.Label.FontWeight='bold';
cb.Title.String='Hz';
cb.Title.FontSize=17;
cb.Title.FontWeight='bold';
local_save(fig,figureRoot,stem);

figPath=figure('Visible','off','Color','w','Position',[80 80 1280 1020]);
heatmapImage=imagesc(betaXValues,betaYValues,plotHC);
heatmapImage.AlphaData=isfinite(plotHC);
set(gca,'YDir','normal','Color',[0.80 0.80 0.80]);
axis xy tight
hold on
plot(betaXValues,bestBetaY,'k-','LineWidth',2.5);
xlabel(xLabel,'FontSize',22,'FontWeight','bold');
ylabel(yLabel,'FontSize',22,'FontWeight','bold');
set(gca,'FontSize',19,'FontWeight','bold','LineWidth',1.25, ...
    'XGrid','off','YGrid','off');
colormap(figPath,jet);
ax=gca;
ax.ColorScale='log';
if finiteValues(lowerIndex)<finiteValues(upperIndex)
    clim([finiteValues(lowerIndex) finiteValues(upperIndex)]);
end
cb=colorbar;
cb.FontSize=17;
cb.Label.String='HC norm';
cb.Label.FontSize=19;
cb.Label.FontWeight='bold';
cb.Title.String='Hz';
cb.Title.FontSize=17;
cb.Title.FontWeight='bold';
local_save(figPath,figureRoot,[stem '_with_Best_Compensation_Path']);

save(fullfile(outputRoot,'full_grid_results.mat'),'gridData','pairType', ...
    'betaXValues','betaYValues','hcMatrix','validHC','plotHC','convergedMatrix', ...
    'residualMatrix','iterationMatrix','bestPath','metrics','-v7.3');
fprintf('Aggregated %s: %d points, %d converged, %d nonconverged.\n', ...
    pairType,totalCount,nnz(convergedMatrix),nnz(~convergedMatrix));
end

function local_save(fig,root,stem)
savefig(fig,fullfile(root,[stem '.fig']));
exportgraphics(fig,fullfile(root,[stem '.pdf']),'ContentType','vector');
close(fig);
end
