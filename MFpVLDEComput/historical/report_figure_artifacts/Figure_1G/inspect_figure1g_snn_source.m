function inspect_figure1g_snn_source(repoRoot)
% Convert the authoritative live script and report its MAT-file contracts.
driverMlx = fullfile(repoRoot, 'Paper2_Fig7Comp_NW_LDE.mlx');
driverM = fullfile(repoRoot, 'Paper2_Fig7Comp_NW_LDE.m');
matlab.internal.liveeditor.openAndConvert(driverMlx, driverM);

fprintf('PARAMVARS\n');
whos('-file', fullfile(repoRoot, 'Data', 'Paper2_NetworkTuning', ...
    'Fig1V4', 'AllMFPixPara_Paper2TuneFig1V4D2.mat'));
fprintf('INIVARS\n');
whos('-file', fullfile(repoRoot, 'Data', 'Initials.mat'));
end
