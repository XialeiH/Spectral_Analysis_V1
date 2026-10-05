function prepare_figure1g_cg_response(runRoot)
% Convert the Paper 3 response table once, outside all timed trials.
repoRoot = fullfile(runRoot, 'repo');
sourceFile = fullfile(repoRoot, 'Data', 'Paper2_NetworkTuning', 'Fig1V4', ...
    'Paper3PlotingData', 'Func200V4D2-L6RealCtrlL4-Ang0.0-SF2.5TF10-Small.mat');
D = load(sourceFile, 'L4EmeshX', 'L4ImeshY', 'LDEFrfunc');
l6Mesh = (1:size(D.LDEFrfunc.S,2)).';
lgnCount = size(D.LDEFrfunc.S,1);
func.S = zeros(lgnCount,300,400,numel(l6Mesh));
func.C = zeros(lgnCount,300,400,numel(l6Mesh));
func.I = zeros(lgnCount,300,400,numel(l6Mesh));
for lgn = 1:lgnCount
    for layer = 1:numel(l6Mesh)
        func.S(lgn,:,:,layer) = D.LDEFrfunc.S{lgn,layer};
        func.C(lgn,:,:,layer) = D.LDEFrfunc.C{lgn,layer};
        func.I(lgn,:,:,layer) = D.LDEFrfunc.I{lgn,layer};
    end
end
l4EMesh = D.L4EmeshX;
l4IMesh = D.L4ImeshY;
save(fullfile(runRoot,'bundles','cg_response_c100_a0.mat'), ...
    'l4EMesh','l4IMesh','l6Mesh','func','-v7.3');
end
