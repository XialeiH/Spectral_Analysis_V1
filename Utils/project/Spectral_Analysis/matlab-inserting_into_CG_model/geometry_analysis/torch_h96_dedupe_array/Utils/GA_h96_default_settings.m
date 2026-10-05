function cfg = GA_h96_default_settings(projectRoot)
% Defaults copied from eigenvalue_5D_mlp_contrast_loop_h96.mlx.
if nargin < 1 || isempty(projectRoot)
    projectRoot = pwd;
end

cfg.ProjectRoot = projectRoot;
cfg.RuntimeDir = [[repro_paths('project') '/Spectral_Analysis/'] ...
    'sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/' ...
    'hpc_runtime_bundle'];
cfg.SaveToFolder = [projectRoot '/Data/Paper2_NetworkTuning/Fig1V4/Paper3PlotingData/'];
cfg.DataFolder = [projectRoot '/Data/Paper2_NetworkTuning/Fig1V4/25Function_Binocular/'];
cfg.DataFolder1 = [projectRoot '/Data/Paper2_NetworkTuning/Fig1V4/16Function_Scheme/LARGE16/'];
cfg.DataFolder2 = [projectRoot '/Data/Paper2_NetworkTuning/Fig1V4/25Function_binocular_realLGN/'];
cfg.DataFolder3 = [projectRoot '/Data/Paper2_NetworkTuning/Fig1V4/Paper2PlotingData/Typical_trajs/'];
cfg.DataFolder4 = [projectRoot '/Data/Paper2_NetworkTuning/Fig1V4/Paper3ICTestData/'];
cfg.FigurePaperPath = [projectRoot '/Figures/Demo062324/'];
cfg.DataPt = 'V4D2';
cfg.MonocuFlag = false;

cfg.xEAll = 0:0.2:1;
cfg.xIAll = 0:0.2:1;
cfg.xEInd = 6;
cfg.xIInd = 6;
cfg.lgnSF = 2.5;
cfg.lgnTF = 10;
cfg.ExpTex = 'CtrlL4';
cfg.Comment = ['L6Real' cfg.ExpTex];
cfg.DomList = {'Small','Large'};

cfg.AngleList = {'0.00','3.75','7.50','11.25','15.00','18.75','22.50'};
cfg.AngleTestAll = 0:3.75:22.5;
cfg.ContrList = [19 66 100];
cfg.TestExSeq = cfg.xEAll(cfg.xEInd) * ones(size(cfg.AngleTestAll));
cfg.TestEySeq = ones(size(cfg.TestExSeq));
cfg.TestIxSeq = cfg.xIAll(cfg.xIInd) * ones(size(cfg.AngleTestAll));
cfg.TestIySeq = ones(size(cfg.TestExSeq));

cfg.N_HCOutX = 4;
cfg.N_HCOutY = 4;
cfg.N_HCOut = 4;
cfg.NPixX = 10;
cfg.NPixY = 10;
cfg.PixNumOut = cfg.N_HCOutX * cfg.N_HCOutY * cfg.NPixX * cfg.NPixY;
cfg.p = 0.33;
cfg.EpocPrecompute = 300;
cfg.EpocTest = 100;
cfg.ExportFlag = 'xn';
cfg.Isaturation = true;

L6pars = {};
L6pars{3} = {2.5,3,[28; 40; 50],[62.5; 84; 96],'quadratic'};
L6pars{1} = {2.5,3,62.5,86,86};
L6pars{2} = {2.5,3,62.5,86,80.4,106,90,142.2, 102, 166.7, 108, 210, 118};

L6ParamRawUse = {2.5,3,62.5,86,80.4,106,90,142.2,102,166.7,108,210,118,{[0 15 29],[0,34.0],[28 32 50],[69 85]}};
L6C1Grid = 0:0.25:120;
L6ParamUse = {2.5,3,{'c1smooth',L6ParamRawUse,L6C1Grid}};
L6pars{5} = L6ParamUse;

cfg.L6pars = L6pars;
cfg.L6parId = 5;
cfg.L6ParamRawUse = L6ParamRawUse;
cfg.L6C1Grid = L6C1Grid;
cfg.L6ParamUse = L6ParamUse;
cfg.L6BaselineChoice = 'c1smooth_LowHighQuadr5_y0_34p0';

IKp7.Thrsld1 = 54;
IKp7.Thrsld2 = 63;
IKp7.Highist = 103;
IKp7.Slope = 0.941;
IKp7.down1 = 0.3;
IKp7.IntaH = 4.5;
IKp7.IntaL = -2.5;
IKp7.IntbL = 4;
IKp7.IntbH = -0.3;
IKp7.Mode = 'multisigmoid';

cfg.IKpUse = IKp7;
cfg.IKp7 = IKp7;
cfg.EKpUse = {1,0,[50; 70; 100],[50; 60.5; 69],'quadratic'};

cfg.HCNorm.wC = 0.3077;
cfg.HCNorm.aE = 32/(32+8);
cfg.HCNorm.aI = 8/(32+8);
end
