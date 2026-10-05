function replot_hc5_snapshot_finite_time_transient(dataFile,setupFile,pdfFile)
% Replot the HC=5 instantaneous-versus-finite-time comparison clearly.

loaded = load(dataFile,'result');
result = loaded.result;
setupLoaded = load(setupFile,'setup');
context = setupLoaded.setup.Context;
mapSize = double(context.MapSize(:).');

snapshotInputE = local_e_map(result.SnapshotInput,mapSize,context.CWeight);
finiteInputE = local_e_map(result.FiniteTimeOptimalInput,mapSize,context.CWeight);
snapshotOutputE = local_e_map(result.SnapshotOutput,mapSize,context.CWeight);
finiteOutputE = local_e_map(result.FiniteTimeOptimalOutput,mapSize,context.CWeight);
inputDifferenceE = finiteInputE-snapshotInputE;
outputDifferenceE = finiteOutputE-snapshotOutputE;

summary = result.Summary;
peakTimeMs = summary.peakTimeMs;
peakGain = summary.finiteTimePeakGain;
inputCosine = summary.inputAbsoluteCosine;
outputCosine = summary.outputAbsoluteCosine;
snapshotSigma = summary.snapshotDPhiSigma1;

figureHandle = figure('Visible','off','Color','w', ...
    'Position',[80 80 1800 920]);
layout = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

inputLimit = max(abs([snapshotInputE(:);finiteInputE(:)]));
outputLimit = max(abs([snapshotOutputE(:);finiteOutputE(:)]));
inputDiffLimit = max(abs(inputDifferenceE(:)));
outputDiffLimit = max(abs(outputDifferenceE(:)));

local_map(nexttile(layout,1),snapshotInputE,inputLimit, ...
    'Current input: top right singular vector of Dphi(x0)');
local_map(nexttile(layout,2),finiteInputE,inputLimit, ...
    sprintf('Finite-time optimal input at Tpeak = %.1f ms',peakTimeMs));
local_map(nexttile(layout,3),inputDifferenceE,inputDiffLimit, ...
    sprintf('Input difference; absolute cosine = %.3f',inputCosine));
axisGain = nexttile(layout,4);
plot(axisGain,result.HorizonGridMs,result.FiniteTimeGain,'-o', ...
    'LineWidth',1.8,'MarkerFaceColor',[0.00 0.45 0.74]);
hold(axisGain,'on');
xline(axisGain,peakTimeMs,'--r',sprintf('Tpeak = %.1f ms',peakTimeMs), ...
    'LabelVerticalAlignment','bottom','Interpreter','none');
xlabel(axisGain,'Finite-time horizon T (ms)');
ylabel(axisGain,'Maximum finite-time gain');
title(axisGain,'Trajectory-aware tangent amplification', ...
    'Interpreter','none');
grid(axisGain,'on'); box(axisGain,'on');

local_map(nexttile(layout,5),snapshotOutputE,outputLimit, ...
    sprintf('Current output: top left singular vector; sigma1(Dphi) = %.2f',snapshotSigma));
local_map(nexttile(layout,6),finiteOutputE,outputLimit, ...
    sprintf('Finite-time response at %.1f ms; gain = %.2f',peakTimeMs,peakGain));
local_map(nexttile(layout,7),outputDifferenceE,outputDiffLimit, ...
    sprintf('Output difference; absolute cosine = %.3f',outputCosine));
axisSummary = nexttile(layout,8);
axis(axisSummary,'off');
text(axisSummary,0.02,0.94,{ ...
    'Comparison definition', ...
    sprintf('Base trajectory: random-positive HC norm = %.1f',summary.initialHCnorm), ...
    'Current method: SVD of instantaneous Dphi(x0)', ...
    'Proposed method: SVD of time-ordered P(T,0)', ...
    sprintf('Peak horizon: %.1f ms',peakTimeMs), ...
    sprintf('Peak finite-time gain: %.4g',peakGain), ...
    sprintf('Input full-state absolute cosine: %.4f',inputCosine), ...
    sprintf('Output full-state absolute cosine: %.4f',outputCosine), ...
    sprintf('Tangent time step: %.1f ms',summary.tangentStepMs)}, ...
    'Units','normalized','VerticalAlignment','top','FontSize',11, ...
    'Interpreter','none');

title(layout,{ ...
    'HC=5 transient direction: instantaneous Jacobian SVD versus finite-time propagator SVD', ...
    'E population maps; contrast 100, orientation 0.00 deg, baseline dynamic L6'}, ...
    'FontWeight','bold','FontSize',14,'Interpreter','none');
exportgraphics(figureHandle,pdfFile,'ContentType','vector');
close(figureHandle);
end

function map = local_e_map(vector,mapSize,wC)
n = prod(mapSize);
map = reshape((1-wC)*real(vector(1:n))+ ...
    wC*real(vector(n+(1:n))),mapSize);
end

function local_map(axisHandle,map,limitValue,titleText)
imagesc(axisHandle,map);
axis(axisHandle,'image');
set(axisHandle,'YDir','normal','XTick',[],'YTick',[]);
colormap(axisHandle,jet(256));
limitValue = max(limitValue,eps);
clim(axisHandle,[-limitValue limitValue]);
colorbar(axisHandle);
title(axisHandle,titleText,'FontSize',10,'Interpreter','none');
end
