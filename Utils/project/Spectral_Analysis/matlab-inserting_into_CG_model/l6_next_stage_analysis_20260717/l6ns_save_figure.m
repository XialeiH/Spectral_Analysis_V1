function l6ns_save_figure(fig, outputDir, stem)
% Save one MATLAB figure in editable, vector, and raster formats.

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
savefig(fig, fullfile(outputDir, [stem '.fig']));
exportgraphics(fig, fullfile(outputDir, [stem '.pdf']), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outputDir, [stem '.png']), 'Resolution', 220);
end
