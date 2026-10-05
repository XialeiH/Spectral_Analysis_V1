function aggregate_real_tuning_results()
% Aggregate beta sweeps and produce requested MATLAB figures.

outputRoot = getenv('REAL_TUNING_OUTPUT_ROOT');
if isempty(outputRoot) || ~exist(outputRoot,'dir')
    error('RealTuning:OutputRoot','REAL_TUNING_OUTPUT_ROOT is required.');
end
figureRoot = fullfile(outputRoot,'figures');
if ~exist(figureRoot,'dir'); mkdir(figureRoot); end
setupFile = getenv('REAL_TUNING_SETUP_FILE');
loadedSetup = load(setupFile,'setup');
contrastValues = unique(loadedSetup.setup.Context.ContrastUse(:));
orientationValues = unique(loadedSetup.setup.Context.OrientationUse(:));
conditionText = sprintf('angle %.2f deg, contrast %.4g', ...
    orientationValues(1),contrastValues(1));
files = dir(fullfile(outputRoot,'*','case_summary.tsv'));
if isempty(files)
    error('RealTuning:NoCases','No case summaries were found.');
end
tables = cell(numel(files),1);
for index = 1:numel(files)
    tables{index} = readtable(fullfile(files(index).folder,files(index).name), ...
        'FileType','text','Delimiter','\t','TextType','string');
end
summary = sortrows(vertcat(tables{:}),'taskId');
writetable(summary,fullfile(outputRoot,'all_case_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

local_plot_max_real(summary,figureRoot,conditionText);
local_plot_fixed_point_shift(summary,figureRoot,conditionText);

representative = summary(summary.isRepresentative & summary.fixedPointConverged,:);
if isempty(representative)
    manifest = cell2table({ ...
        'maxReal_vs_beta','both',NaN,fullfile(figureRoot, ...
            '01_MaxReal_Eigenvalue_vs_Beta_Real_and_FPP.pdf'); ...
        'fixedPoint_shift','both',NaN,fullfile(figureRoot, ...
            '02_Moved_Fixed_Point_vs_Beta.pdf')}, ...
        'VariableNames',{'figureType','pathway','beta','pdfFile'});
    writetable(manifest,fullfile(outputRoot,'figure_manifest.tsv'), ...
        'FileType','text','Delimiter','\t');
    save(fullfile(outputRoot,'real_tuning_summary.mat'),'summary','representative', ...
        'conditionText','-v7.3');
    fprintf('Aggregated %d smoke cases; no full spectra were requested.\n',height(summary));
    return
end
spectra = cell(height(representative),1);
for index = 1:height(representative)
    match = dir(fullfile(outputRoot,sprintf('%03d_*',representative.taskId(index)), ...
        'case_result.mat'));
    if numel(match)~=1
        error('RealTuning:RepresentativeFile', ...
            'Expected one result file for task %d.',representative.taskId(index));
    end
    resultFile = fullfile(match.folder,match.name);
    loaded = load(resultFile,'eigenvaluesTrue');
    spectra{index} = loaded.eigenvaluesTrue;
end
binWidth = 0.05;
histogramEdges = cell(height(representative),1);
spectrumRanges = nan(height(representative),2);
imaginaryLimits = nan(height(representative),1);

manifestRows = cell(2+2*height(representative),4);
manifestRows(1,:) = {'maxReal_vs_beta','both',NaN, ...
    fullfile(figureRoot,'01_MaxReal_Eigenvalue_vs_Beta_Real_and_FPP.pdf')};
manifestRows(2,:) = {'fixedPoint_shift','both',NaN, ...
    fullfile(figureRoot,'02_Moved_Fixed_Point_vs_Beta.pdf')};
row = 2;
for index = 1:height(representative)
    one = representative(index,:);
    eigenvalues = spectra{index};
    [edges,realRange,imaginaryLimit] = local_spectrum_limits(eigenvalues,binWidth);
    histogramEdges{index} = edges;
    spectrumRanges(index,:) = realRange;
    imaginaryLimits(index) = imaginaryLimit;
    tag = sprintf('%s_Beta_%s',char(one.pathway),local_number_tag(one.beta));
    histogramFile = fullfile(figureRoot,['Histogram_' tag '.pdf']);
    spectrumFile = fullfile(figureRoot,['Eigenspectrum_' tag '.pdf']);
    local_plot_histogram(eigenvalues,edges,one,conditionText,histogramFile);
    local_plot_spectrum(eigenvalues,realRange,imaginaryLimit, ...
        one,conditionText,spectrumFile);
    row = row+1;
    manifestRows(row,:) = {'histogram',char(one.pathway),one.beta,histogramFile};
    row = row+1;
    manifestRows(row,:) = {'eigenspectrum',char(one.pathway),one.beta,spectrumFile};
end
manifest = cell2table(manifestRows,'VariableNames', ...
    {'figureType','pathway','beta','pdfFile'});
writetable(manifest,fullfile(outputRoot,'figure_manifest.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputRoot,'real_tuning_summary.mat'),'summary','representative', ...
    'spectra','histogramEdges','spectrumRanges','imaginaryLimits','binWidth', ...
    'conditionText','-v7.3');
fprintf('Aggregated %d cases and generated %d figures.\n',height(summary),height(manifest));
end

function local_plot_max_real(summary,figureRoot,conditionText)
fig = figure('Visible','off','Color','w','Position',[100 100 1450 590]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
names = ["L6","Inhibition"];
for index = 1:2
    ax = nexttile(layout,index);
    one = sortrows(summary(summary.pathway==names(index),:),'beta');
    plot(ax,one.beta,one.trueMaxRealLambda,'-o','LineWidth',2, ...
        'MarkerSize',4,'DisplayName','real tuning, moved fixed point'); hold(ax,'on');
    plot(ax,one.beta,one.fppMaxRealLambda,'--s','LineWidth',1.8, ...
        'MarkerSize',4,'DisplayName','FPP, matched w=-beta');
    yline(ax,1,'k:','LineWidth',1.5,'DisplayName','stability boundary');
    grid(ax,'on'); xlabel(ax,'\beta (fractional whole-pathway change)');
    ylabel(ax,'max Re \lambda(D_x\Phi)');
    title(ax,sprintf('%s pathway',names(index)),'FontWeight','bold');
    legend(ax,'Location','best');
end
title(layout,sprintf(['Regular iteration after real whole-pathway tuning versus ' ...
    'fixed-point preserving tuning: %s'],conditionText), ...
    'FontWeight','bold');
local_save(fig,figureRoot,'01_MaxReal_Eigenvalue_vs_Beta_Real_and_FPP');
end

function local_plot_fixed_point_shift(summary,figureRoot,conditionText)
fig = figure('Visible','off','Color','w','Position',[100 100 1450 590]);
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
names = ["L6","Inhibition"];
for index = 1:2
    ax = nexttile(layout,index);
    one = sortrows(summary(summary.pathway==names(index),:),'beta');
    yyaxis(ax,'left');
    plot(ax,one.beta,one.relativeFixedPointShift,'-o','LineWidth',2,'MarkerSize',4);
    ylabel(ax,'||f^*_{\beta}-f^*_0||_2 / ||f^*_0||_2');
    yyaxis(ax,'right');
    plot(ax,one.beta,one.HCnormFromBaseline,'-s','LineWidth',1.8,'MarkerSize',4);
    ylabel(ax,'HCnorm from baseline');
    grid(ax,'on'); xlabel(ax,'\beta (fractional whole-pathway change)');
    title(ax,sprintf('%s pathway',names(index)),'FontWeight','bold');
end
title(layout,sprintf('Operating-point displacement produced by real pathway tuning: %s', ...
    conditionText), ...
    'FontWeight','bold');
local_save(fig,figureRoot,'02_Moved_Fixed_Point_vs_Beta');
end

function local_plot_histogram(eigenvalues,edges,one,conditionText,outputFile)
fig = figure('Visible','off','Color','w','Position',[100 100 1050 650]);
histogram(real(eigenvalues),edges,'FaceColor',[0.18 0.52 0.74], ...
    'EdgeColor','none'); grid on; xlim([edges(1) edges(end)]);
if edges(1)<=1 && edges(end)>=1; xline(1,'r--','LineWidth',1.5); end
xlabel('Re \lambda(D_x\Phi)'); ylabel('eigenvalue count per 0.05 bin');
title(sprintf(['%s real tuning: beta=%+.3f, gain=%.3f, moved max Re lambda=%.6f\n' ...
    '%s; fixed-point shift=%.4g, residual=%.2e, real-axis range=[%.2f, %.2f]'], ...
    one.pathway,one.beta,one.pathwayGain,one.trueMaxRealLambda, ...
    conditionText,one.relativeFixedPointShift,one.fixedPointResidual, ...
    edges(1),edges(end)), ...
    'FontWeight','bold');
exportgraphics(fig,outputFile,'ContentType','vector'); close(fig);
end

function local_plot_spectrum(eigenvalues,realRange,imaginaryLimit,one,conditionText,outputFile)
fig = figure('Visible','off','Color','w','Position',[100 100 1050 650]);
scatter(real(eigenvalues),imag(eigenvalues),12,[0.10 0.42 0.70],'filled', ...
    'MarkerFaceAlpha',0.60); grid on; hold on;
if realRange(1)<=1 && realRange(2)>=1; xline(1,'r--','LineWidth',1.5); end
yline(0,'k:');
xlim(realRange); ylim([-imaginaryLimit imaginaryLimit]);
xlabel('Re \lambda(D_x\Phi)'); ylabel('Im \lambda(D_x\Phi)');
title(sprintf(['%s real tuning eigenspectrum: beta=%+.3f, gain=%.3f, ' ...
    'moved max Re lambda=%.6f\n%s; real-axis range matches its 0.05-bin histogram: [%.2f, %.2f]'], ...
    one.pathway,one.beta,one.pathwayGain,one.trueMaxRealLambda, ...
    conditionText,realRange(1),realRange(2)),'FontWeight','bold');
exportgraphics(fig,outputFile,'ContentType','vector'); close(fig);
end

function [edges,realRange,imaginaryLimit] = local_spectrum_limits(eigenvalues,binWidth)
realMinimum = floor(min(real(eigenvalues))/binWidth)*binWidth;
realMaximum = ceil(max(real(eigenvalues))/binWidth)*binWidth;
if realMaximum<=realMinimum; realMaximum=realMinimum+binWidth; end
if max(real(eigenvalues))>0.5
    realMaximum = max(realMaximum,1.05);
end
edges = realMinimum:binWidth:realMaximum;
if edges(end)<realMaximum; edges(end+1)=realMaximum; end
realRange = [edges(1),edges(end)];
imaginaryLimit = max(1e-4,1.05*max(abs(imag(eigenvalues))));
end

function local_save(fig,root,stem)
savefig(fig,fullfile(root,[stem '.fig']));
exportgraphics(fig,fullfile(root,[stem '.pdf']),'ContentType','vector');
close(fig);
end

function tag = local_number_tag(value)
tag = sprintf('%+.3f',value);
tag = strrep(tag,'+','p'); tag = strrep(tag,'-','m'); tag = strrep(tag,'.','p');
end
