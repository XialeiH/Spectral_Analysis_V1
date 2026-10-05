function make_Final_Figure_3E_2()
% Export the three analytical Figure 3F panels separately.

projectRoot = '/Users/xialeihuang/Desktop/Neuroscience_Project';
modelRoot = fullfile(projectRoot, 'Spectral_Analysis', ...
    'matlab-inserting_into_CG_model');
setupFile = fullfile(modelRoot, 'l6_next_stage_analysis_20260717', ...
    'results_global_bifurcation_20260722', ...
    'global_bifurcation_setup.mat');
outputRoot = tempdir;
outputFig = fullfile(projectRoot, 'Spectral_Analysis_Report', 'Figures', ...
    'Final draft', 'Figure 3E.2.fig');

loaded = load(setupFile, 'setup');
J = loaded.setup.Pathway.JBaseline;
n = size(J, 1) / 3;
e = 1:(2*n);
i = (2*n+1):(3*n);
A = J(e,e);
B = J(e,i);
C = J(i,e);
D = J(i,i);
returnedPath = C * (A \ B);
H = D - returnedPath;

fprintf('Computing the 1,600-state Schur eigensystem.\n');
[V, lambda] = eig(full(H), 'vector');
V = V ./ max(vecnorm(V, 2, 1), eps);
directAction = D * V;
returnedAction = returnedPath * V;
a = vecnorm(directAction, 2, 1).';
b = vecnorm(returnedAction, 2, 1).';
r = vecnorm(directAction - returnedAction, 2, 1).' ./ ...
    max(a + b, eps);
kappa = spatial_frequency_score(V, [40 40]);

neuroColors = [
    178  24  43
    227  74  51
    254 227 145
    171 217 233
     67 147 195
     33 102 172] / 255;
deepRed = neuroColors(1,:);
deepBlue = neuroColors(6,:);
guideGrey = [0.55 0.55 0.55];
continuousColors = interp1(linspace(0,1,size(neuroColors,1)), ...
    neuroColors, linspace(0,1,256), 'linear');

% Figure 3F.2: pathway balance.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 10.5 8.5]);
ax = axes(fig, 'Position', [0.16 0.16 0.66 0.72]);
hold(ax, 'on');
scatter(ax, a, b, 32, r, 'filled', ...
    'MarkerFaceAlpha', 0.80, 'MarkerEdgeAlpha', 0.80);
positiveAB = [a(a>0); b(b>0)];
abLimits = [max(min(positiveAB)*0.8, 1e-8), max([a;b])*1.25];
plot(ax, abLimits, abLimits, '--', 'Color', guideGrey, ...
    'LineWidth', 2.5, 'HandleVisibility', 'off');
set(ax, 'XScale', 'log', 'YScale', 'log');
xlim(ax, abLimits); ylim(ax, abLimits);
axis(ax, 'square');
xlabel(ax, '$a_k=\|Dv_k\|_2$', 'Interpreter', 'latex');
ylabel(ax, '$b_k=\|CA^{-1}Bv_k\|_2$', 'Interpreter', 'latex');
title(ax, 'Mode-resolved pathway balance', ...
    'FontSize', 29, 'FontWeight', 'bold');
colormap(ax, continuousColors);
clim(ax, [0 1]);
cb = colorbar(ax, 'eastoutside');
cb.Position = [0.845 0.16 0.030 0.72];
cb.Title.String = '$r_k$';
cb.Title.Interpreter = 'latex';
cb.FontSize = 26;
cb.Ticks = 0:0.2:1;
format_axes(ax);
savefig(fig, outputFig);
exportgraphics(fig, fullfile(outputRoot, 'Figure 3F.2.pdf'), ...
    'ContentType', 'vector');
close(fig);

% Figure 3F.3: weak spectral modes versus cancellation.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 9.5 7.5]);
ax = axes(fig, 'Position', [0.17 0.17 0.76 0.72]);
scatter(ax, max(r, eps), max(abs(lambda), eps), 32, ...
    'MarkerFaceColor', deepBlue, 'MarkerEdgeColor', deepBlue, ...
    'MarkerFaceAlpha', 0.76, 'MarkerEdgeAlpha', 0.76);
set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, '$\mathrm{Cancellation\ residual}\ r_k$', ...
    'Interpreter', 'latex');
ylabel(ax, '$|\lambda_k(H)|$', 'Interpreter', 'latex');
title(ax, 'Weak modes and cancellation', ...
    'FontSize', 29, 'FontWeight', 'bold');
format_axes(ax);
exportgraphics(fig, fullfile(outputRoot, 'Figure 3F.3.pdf'), ...
    'ContentType', 'vector');
close(fig);

% Figure 3F.4: cancellation versus spatial frequency.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 9.5 7.5]);
ax = axes(fig, 'Position', [0.17 0.17 0.76 0.72]);
scatter(ax, kappa, max(r, eps), 32, ...
    'MarkerFaceColor', deepRed, 'MarkerEdgeColor', deepRed, ...
    'MarkerFaceAlpha', 0.76, 'MarkerEdgeAlpha', 0.76);
set(ax, 'YScale', 'log');
xlabel(ax, ...
    '$\mathrm{Spatial\ frequency}\ \kappa_k\ \mathrm{(cycles/pixel)}$', ...
    'Interpreter', 'latex');
ylabel(ax, '$\mathrm{Cancellation\ residual}\ r_k$', ...
    'Interpreter', 'latex');
title(ax, 'Cancellation within spatial modes', ...
    'FontSize', 29, 'FontWeight', 'bold');
format_axes(ax);
exportgraphics(fig, fullfile(outputRoot, 'Figure 3F.4.pdf'), ...
    'ContentType', 'vector');
close(fig);

fprintf('Saved Figure 3F.2.pdf, Figure 3F.3.pdf, and Figure 3F.4.pdf.\n');
end

function score = spatial_frequency_score(V, mapSize)
nModes = size(V, 2);
rows = mapSize(1);
columns = mapSize(2);
fx = (-floor(columns/2):ceil(columns/2)-1) / columns;
fy = (-floor(rows/2):ceil(rows/2)-1) / rows;
[kx, ky] = meshgrid(fx, fy);
radiusSquared = kx.^2 + ky.^2;
score = zeros(nModes, 1);
for modeIndex = 1:nModes
    spatialMap = reshape(V(:,modeIndex), mapSize);
    energy = abs(fftshift(fft2(spatialMap))).^2;
    score(modeIndex) = sqrt(sum(radiusSquared .* energy, 'all') / ...
        max(sum(energy, 'all'), eps));
end
end

function format_axes(ax)
set(ax, 'FontName', 'Arial', 'FontSize', 26, 'FontWeight', 'bold', ...
    'LineWidth', 1.4, 'Box', 'on', 'Layer', 'top', ...
    'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'off');
ax.Toolbar.Visible = 'off';
end
