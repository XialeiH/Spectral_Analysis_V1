tic

ThisScriptDir = repro_paths('project/Spectral_Analysis/matlab-inserting_into_CG_model/l6_next_stage_analysis_20260717');
if isempty(ThisScriptDir)
    ThisScriptDir = pwd;
end
CanonicalCodeDir = getenv('L6NS_CANONICAL_CODE_DIR');
if ~isempty(CanonicalCodeDir)
    ThisScriptDir = CanonicalCodeDir;
end
TaskId = str2double(getenv('SLURM_ARRAY_TASK_ID'));
if isnan(TaskId)
    TaskId = str2double(getenv('GEOM_TASK_ID'));
end
if isnan(TaskId)
    TaskId = 1;
end
AngleGrid = [0:1:22, 22.5];
ContrGrid = 2:2:100;
if TaskId < 1 || TaskId > numel(AngleGrid) * numel(ContrGrid)
    error('TaskId %d outside valid range 1:%d.', TaskId, numel(AngleGrid) * numel(ContrGrid));
end
[AngleGridInd, ContrGridInd] = ind2sub([numel(AngleGrid), numel(ContrGrid)], TaskId);
AngleCond = AngleGrid(AngleGridInd);
ContrCond = ContrGrid(ContrGridInd);
EnvAngleCond = getenv('GEOM_ANGLE_COND');
if ~isempty(EnvAngleCond)
    AngleCond = str2double(EnvAngleCond);
end
EnvContrCond = getenv('GEOM_CONTR_COND');
if ~isempty(EnvContrCond)
    ContrCond = str2double(EnvContrCond);
end

MainRoot = getenv('GEOM_MAIN_ROOT');
if isempty(MainRoot)
    MainRoot = '/scratch/xh2906/NYU-Vision-2Drive-main';
end
GeometryRoot = getenv('GEOM_ROOT');
if isempty(GeometryRoot)
    GeometryRoot = ThisScriptDir;
end
RunRoot = getenv('GEOM_RUN_ROOT');
if isempty(RunRoot)
    RunRoot = fullfile(GeometryRoot, 'results', 'full_h96_geometry_torch_array');
end
CurrentFolder = MainRoot
repro_addpath(CurrentFolder)
repro_addpath([CurrentFolder '/Utils'])
OverrideUtils = fullfile(ThisScriptDir, 'Utils');
if exist(OverrideUtils, 'dir')
    repro_addpath(OverrideUtils, '-begin')
end
h96RuntimeDir = getenv('H96_RUNTIME_DIR');
if isempty(h96RuntimeDir)
    h96RuntimeDir = '/scratch/xh2906/librarySCI_runs/sixd_prefAngle_varaware_singlemlp_h96_final70_smooth2d_low005_stratmid_torch_20260518_044134/hpc_runtime_bundle';
end
repro_addpath(h96RuntimeDir, '-begin')
if exist(OverrideUtils, 'dir')
    repro_addpath(OverrideUtils, '-begin')
end
repro_addpath([CurrentFolder '/Data'])
SaveToFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper3PlotingData/']; % V1D2
DataFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/25Function_Binocular/']; % V1D2
DataFolder1 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/16Function_Scheme/LARGE16/']; % V1D2
DataFolder2 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/25Function_binocular_realLGN/']; % V1D2
DataFolder3 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper2PlotingData/Typical_trajs/']; % V1D2
DataFolder4 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper3ICTestData/'];
repro_addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'])

GeometryOutputRoot = fullfile(RunRoot, 'geometry_mat');
FigureOutputRoot = fullfile(RunRoot, 'figures', 'spectral_analysis_eigenvalue_eigenvectors');
if ~exist(GeometryOutputRoot, 'dir')
    mkdir(GeometryOutputRoot)
end
if ~exist(FigureOutputRoot, 'dir')
    mkdir(FigureOutputRoot)
end

repro_addpath(DataFolder4)
repro_addpath(DataFolder3)
repro_addpath(SaveToFolder)
repro_addpath(DataFolder2)
repro_addpath(DataFolder1)
repro_addpath(DataFolder)
repro_addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'])
fprintf('Torch array task %d: angle %.2f, contrast %d\n', TaskId, AngleCond, ContrCond);

DataPt = 'V4D2'; %'V4D1'; 'V5D1'
load(sprintf("AllMFPixPara_Paper2TuneFig1%s.mat",DataPt))

CurrentFolder = MainRoot;
FigurePaperPath = [CurrentFolder '/Figures/Demo062324/']; % V1D2
MonocuFlag = false;
close all

close all
xEAll = 0:0.2:1; xIAll = 0:0.2:1;
lgnSF = 2.5; lgnTF = 10;

ExpTex = 'CtrlL4';
Comment = ['L6Real' ExpTex];

L6pars = {};
L6ParamRawUse = {2.5,3,62.5,86,80.4,106,90,142.2,102,166.7,108,210,118,{[0 15 29],[0,34.0],[28 32 50],[69 85]}};
L6C1Grid = 0:0.25:120;
L6ParamUse = {2.5,3,{'c1smooth',L6ParamRawUse,L6C1Grid}}; % baseline L6 used in iteration
L6AmplitudeScale = 1;
EnvL6AmplitudeScale = getenv('L6_AMPLITUDE_SCALE');
if ~isempty(EnvL6AmplitudeScale)
    L6AmplitudeScale = str2double(EnvL6AmplitudeScale);
    if ~isfinite(L6AmplitudeScale)
        error('L6_AMPLITUDE_SCALE must be numeric.');
    end
end
if L6AmplitudeScale ~= 1
    L6ParamUse = {2.5,3,{'scale',L6AmplitudeScale,L6ParamUse}};
end
L6pars{5} = L6ParamUse;
L6EquBlendWeight = 0; % 0=original dynamic L6; 1=fixed equilibrium L6
EnvL6EquBlendWeight = getenv('L6_EQU_BLEND_WEIGHT');
if ~isempty(EnvL6EquBlendWeight)
    if contains(EnvL6EquBlendWeight, ',')
        L6EquBlendWeight = sscanf(EnvL6EquBlendWeight, '%f,')';
    else
        L6EquBlendWeight = str2double(EnvL6EquBlendWeight);
    end
end
if numel(L6EquBlendWeight) ~= 1 && numel(L6EquBlendWeight) ~= 3
    error('L6EquBlendWeight must be scalar w or vector [wS,wC,wI].');
end
L6TonicScale = 1;
EnvL6TonicScale = getenv('L6_TONIC_SCALE');
L6TonicScaleActive = ~isempty(EnvL6TonicScale);
if L6TonicScaleActive
    L6TonicScale = str2double(EnvL6TonicScale);
    if ~isfinite(L6TonicScale) || L6TonicScale < 0
        error('L6_TONIC_SCALE must be nonnegative numeric.');
    end
end




IKp7.Thrsld1 = 54;
IKp7.Thrsld2 = 63;
IKp7.Highist = 103;
IKp7.Slope   = 0.941;
IKp7.down1   = 0.3;
IKp7.IntaH   = 4.5;
IKp7.IntaL   = -2.5;
IKp7.IntbL   = 4;
IKp7.IntbH   = -0.3;
IKp7.Mode    = 'multisigmoid';

IKpUse = IKp7;

close all

EKpUse = {1,0,[50; 70; 100],[50; 60.5; 69],'quadratic'}; 
Isaturation = true;

close all
%Comment = '';
xEInd = 6%:length(xEAll)
xIInd = 6
fprintf('xE=%.1f, xI=%.1f\n',xEAll(xEInd),xIAll(xIInd))
L6parId = 5
fprintf('L6 ParameterID %d\n',L6parId)
if iscell(L6pars{L6parId}{end}) && numel(L6pars{L6parId}{end}) >= 1 && strcmp(L6pars{L6parId}{end}{1}, 'weightedblend')
    L6Comment = sprintf('L6wbR%.2f', L6pars{L6parId}{end}{2})
elseif ischar(L6pars{L6parId}{end}) && strcmp(L6pars{L6parId}{end}, 'quadratic')
    L6Comment = sprintf('%s',L6pars{L6parId}{end})
elseif iscell(L6pars{L6parId}{end})
    L6Comment = sprintf('LowHighQuadr5')
else
    L6Comment = sprintf('linear')
end
if numel(L6EquBlendWeight) == 1
    L6EquComment = sprintf('L6eqW%.2f', L6EquBlendWeight);
else
    L6EquComment = sprintf('L6eqWS%.2f_C%.2f_I%.2f', L6EquBlendWeight(1), L6EquBlendWeight(2), L6EquBlendWeight(3));
end
if exist('L6AmplitudeScale', 'var') && L6AmplitudeScale ~= 1
    L6EquComment = sprintf('%s_L6amp%.2f', L6EquComment, L6AmplitudeScale);
end
if exist('L6TonicScaleActive', 'var') && L6TonicScaleActive
    L6EquComment = sprintf('%s_L6ton%.2f', L6EquComment, L6TonicScale);
end
if exist('L6KernelSigmaScale', 'var') && L6KernelSigmaScale ~= 1
    L6EquComment = sprintf('%s_L6sig%.2f', L6EquComment, L6KernelSigmaScale);
end
L6EquComment = strrep(L6EquComment, '.', 'p');
L6EquComment = strrep(L6EquComment, '-', 'm');
tic
%AngleList = {'0.00'}; % only vertical for now
AngleList = {sprintf('%.2f', AngleCond)};
DomList =  {'Small','Large'};

ContrList = ContrCond;

% there are 3200 pixels in total
N_HCOutX = 4; N_HCOutY = 4; PixNumOut = N_HCOutX*N_HCOutY*NPixX*NPixY;
N_HCinX = 4; N_HCinY = 4;
NPixX = 10; NPixY = 10;

% here I test multiple angles, 
% For the Sigmoid version, only test 0 degree 
AngleTestAll = AngleCond; % Torch array condition; h96 baseline is continuous in condition space
TestNum = length(AngleTestAll); ICTestAll = cell(TestNum,1);
TestExSeq = xEAll(xEInd)*ones(size(AngleTestAll)); TestEySeq = 1*ones(size(TestExSeq));
TestIxSeq = xIAll(xIInd)*ones(size(AngleTestAll)); TestIySeq = 1*ones(size(TestExSeq));

ICData = load(sprintf('LDETracesV4D2_NewSmear_ang0.0.mat'));
load('LDETracesV4D2_NewSmear_ang0.0.mat')
LDEICPart = ICData.LDEEquv;

LDEICOut = HCRotXY(LDEICPart,0,N_HCinX,N_HCinY,N_HCOutX,N_HCOutY,NPixX,NPixY,0);
for TestId = 1:TestNum
    ICTestAll{TestId} = LDEICOut;
end
LDEOutAll = cell(size(ICTestAll));
NANFlag = false(size(ICTestAll));
%BGFlagAllTest = cell(size(ICTestAll));

p = 0.33;
EpocTest = 100;
ExportFlag = 'xn';

L2Diff_OneStepAll = zeros(length(ICTestAll),EpocTest+1);
DiffVecAll = zeros(length(ICTestAll),EpocTest+1,3*NPixY*NPixX);

% Determine L4 Input from a range
L4SE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4SI = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4CE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4CI = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4IE = zeros(N_HC*NPixX*N_HC*NPixY,1);
L4II = zeros(N_HC*NPixX*N_HC*NPixY,1);
for PInd = 1:N_HC*NPixX*N_HC*NPixY
    L4SE(PInd) = (C_SS_Pixel_Us(PInd,:)*FrSPixVec          + C_SC_Pixel_Us(PInd,:)*FrCPixVec);
    L4SI(PInd) =  C_SI_Pixel_Us(PInd,:)*FrIPixVec ;
    L4CE(PInd) = (C_CS_Pixel_Us(PInd,:)*FrSPixVec          + C_CC_Pixel_Us(PInd,:)*FrCPixVec);
    L4CI(PInd) =  C_CI_Pixel_Us(PInd,:)*FrIPixVec ;
    L4IE(PInd) = (C_IS_Pixel_Us(PInd,:)*FrSPixVec          + C_IC_Pixel_Us(PInd,:)*FrCPixVec);
    L4II(PInd) =  C_II_Pixel_Us(PInd,:)*FrIPixVec ;
end
L4Eall = L4SE+L4CE+L4IE; L4Iall = L4SI+L4CI+L4II;
L4SEp = mean(L4SE./L4Eall);L4CEp = mean(L4CE./L4Eall);L4IEp = mean(L4IE./L4Eall);
L4SIp = mean(L4SI./L4Iall);L4CIp = mean(L4CI./L4Iall);L4IIp = mean(L4II./L4Iall);

for Contr = ContrList
    if Contr == 100
        ContrStr = '';
    else
        ContrStr = sprintf('_Contr%d',Contr);
    end
    fprintf('Running contrast %d\n', Contr);

    EpocTest = 300;
    LDEOutAll = cell(size(ICTestAll));
    LDEOut = cell(size(ICTestAll));
    NANFlag = false(size(ICTestAll));
    L2Diff_OneStepAll = zeros(length(ICTestAll),EpocTest+1);
    DiffVecAll = zeros(length(ICTestAll),EpocTest+1,3*NPixY*NPixX);
    err_ode_true_all = zeros(length(ICTestAll),EpocTest+1);

% test every h96 baseline condition
for TestId = 1:TestNum
AngleInpt = AngleTestAll(TestId);

RotInd = floor(AngleInpt/45); % 0-3
MirInd = mod(floor(AngleInpt/22.5),2); % 0: no mirror; 1: mirror
switch MirInd % find the corresponding source
    case 0
        AngSource = mod(AngleInpt,45);
    case 1
        AngSource = mod(-AngleInpt,45);
end

% h96 baseline accepts continuous alpha/contrast; keep the exact array angle.
AngFuncCtgrCdid = [0:1:22, 22.5];
[~,AngFuncCtgr] = min(abs(AngFuncCtgrCdid - AngleInpt));
AngPrint = sprintf('%.2f', AngleInpt);
Angle = AngleInpt;

L4EmeshXAll = {[]};
L4ImeshYAll = {[]};
LDEFrfuncAll = {struct('S',{{}},'C',{{}},'I',{{}})};

%    N_HCOut = 4; PixNumOut = N_HCOut^2*NPixX*NPixY;
LGNlist = 5; L6list = 33; % Why do we have L6list here?
OD_SMapModi = OD_SMap;


%OD_SMapModi = OD_SMap(1:2*n_S_HC,1:2*n_S_HC);
if MonocuFlag
    OD_SMapModi(n_S_HC+1:2*n_S_HC,1:end) = 5;
    SaveStr = ['Mono' Comment];
else
    SaveStr = ['Bino' Comment];
end

L6ConvSig = 0.75;
L6KernelSigmaScale = 1;
EnvL6KernelSigmaScale = getenv('L6_KERNEL_SIGMA_SCALE');
if ~isempty(EnvL6KernelSigmaScale)
    L6KernelSigmaScale = str2double(EnvL6KernelSigmaScale);
    if ~isfinite(L6KernelSigmaScale) || L6KernelSigmaScale <= 0
        error('L6_KERNEL_SIGMA_SCALE must be positive numeric.');
    end
end
if numel(L6EquBlendWeight) == 1
    L6EquComment = sprintf('L6eqW%.2f', L6EquBlendWeight);
else
    L6EquComment = sprintf('L6eqWS%.2f_C%.2f_I%.2f', L6EquBlendWeight(1), L6EquBlendWeight(2), L6EquBlendWeight(3));
end
if exist('L6AmplitudeScale', 'var') && L6AmplitudeScale ~= 1
    L6EquComment = sprintf('%s_L6amp%.2f', L6EquComment, L6AmplitudeScale);
end
if exist('L6TonicScaleActive', 'var') && L6TonicScaleActive
    L6EquComment = sprintf('%s_L6ton%.2f', L6EquComment, L6TonicScale);
end
if exist('L6KernelSigmaScale', 'var') && L6KernelSigmaScale ~= 1
    L6EquComment = sprintf('%s_L6sig%.2f', L6EquComment, L6KernelSigmaScale);
end
L6EquComment = strrep(L6EquComment, '.', 'p');
L6EquComment = strrep(L6EquComment, '-', 'm');

L6ConvSig = L6ConvSig * L6KernelSigmaScale;
L6ConvTruc = 1.5/(L6ConvSig);
Kersize = floor(2*L6ConvSig*L6ConvTruc);
if Kersize<=0
    Kersize= 2; % at least 1
end
[x,y] = meshgrid(1:Kersize,1:Kersize);
c = Kersize/2;
exponent = ((x-c).^2+(y-c).^2)/(2*L6ConvSig^2);
L6Kernel = exp(-exponent);
L6Kernel(L6Kernel<exp(-(L6ConvTruc*1.0)^2/2))=0; % truncate
L6Kernel = L6Kernel + ...
    L6Kernel(end:-1:1,:) + ...
    L6Kernel(:,end:-1:1) + ...
    L6Kernel(end:-1:1,end:-1:1);%symmetrize
L6Kernel = L6Kernel/sum(L6Kernel,'all'); % normalize
L6Kernel = L6Kernel/((sum(L6Kernel,'all') - L6Kernel(2,2))/(1-0.3));
L6Kernel(2,2) = 0.3;

N_HCs = 3; % shrink to 2*2
FigOn = false;
LGNsmear = n_S_HC*0.2;
TruncLGN = 1;
LGNFilt_Grating = SpatialGaussianFilt_my(OD_SMapModi,...
    N_HCs,n_S_HC,LGNsmear,TruncLGN,FigOn);

if ~MonocuFlag
    LGNFilt_Grating = [LGNFilt_Grating,zeros(size(LGNFilt_Grating,1),1)];
end

PixLGNCtgr = LGNIndSpat_Rec(...
    LGNFilt_Grating,1:LGNlist,NnSPixel,N_HC,N_HCOutX,N_HCOutY,NPixX,NPixY,true);

% symmetrize !! This is actually trick, so I do it ad hoc
PixLGNCtgr = Ocu_LGNL6symm(PixLGNCtgr,N_HCOutX,N_HCOutY,NPixX,NPixY);

BGFlag = false;
if BGFlag
    PixLGNCtgr(:,end) = 1; PixLGNCtgr(:,1:end-1) = 0;
end

MirInd = logical(MirInd);
CtgrOrder = [2,1;3,4]; % for lower left HC
CtgrOrderUse = rot90(CtgrOrder,RotInd);
if MirInd
    %flipDim = mod(RotInd,2)+1; %1 for col 2 ofr row
    CtgrOrderUse = flip(CtgrOrderUse,1);
end
CtgrOrderReadout = ...
    [CtgrOrderUse(1,2);
    CtgrOrderUse(1,1);
    CtgrOrderUse(2,1);
    CtgrOrderUse(2,2)];
%    for Ctgr = 1:4
%         PixL6Ctgr(:,Ctgr) = ...
%             HCRot_Rec(PixL6Ctgr(:,Ctgr),RotInd,N_HCOutX,N_HCOutY,NPixX,NPixY,MirInd);
%         PixLGNCtgr(:,Ctgr) = ...
%             HCRot_Rec(PixLGNCtgr(:,Ctgr),RotInd,N_HCOutX,N_HCOutY,NPixX,NPixY,MirInd);
%    end

% let's implement something cheap now for 45 deg
PixLGNCtgr = PixLGNCtgr(:,[CtgrOrderReadout',5]);

ParaEx = TestExSeq(TestId); ParaEy = TestEySeq(TestId);
ParaIx = TestIxSeq(TestId); ParaIy = TestIySeq(TestId);

C_SS_mean = sparse(AveSpatKer_Rec(C_SS_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
C_CS_mean = sparse(AveSpatKer_Rec(C_CS_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
C_IS_mean = sparse(AveSpatKer_Rec(C_IS_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
C_SC_mean = sparse(AveSpatKer_Rec(C_SC_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
C_CC_mean = sparse(AveSpatKer_Rec(C_CC_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
C_IC_mean = sparse(AveSpatKer_Rec(C_IC_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaEx,ParaEy));
C_SI_mean = sparse(AveSpatKer_Rec(C_SI_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaIx,ParaIy));
C_CI_mean = sparse(AveSpatKer_Rec(C_CI_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaIx,ParaIy));
C_II_mean = sparse(AveSpatKer_Rec(C_II_Pixel_Us,N_HC,...
    N_HCOutX,N_HCOutY,NPixX,NPixY,ParaIx,ParaIy));

Conn_Reduc_TestS = 6; Conn_Reduc_TestC = 4;
Reduc_S = Conn_Reduc_TestS*(1-CplxR); Reduc_C = Conn_Reduc_TestC*CplxR;
P_SS_rate = (N_SS-Reduc_S)/N_SS; C_SS_meanU = C_SS_mean * P_SS_rate;
P_SC_rate = (N_SC-Reduc_C)/N_SC; C_SC_meanU = C_SC_mean * P_SC_rate;
P_CS_rate = (N_CS-Reduc_S)/N_CS; C_CS_meanU = C_CS_mean * P_CS_rate;
P_CC_rate = (N_CC-Reduc_C)/N_CC; C_CC_meanU = C_CC_mean * P_CC_rate;

tic
ICUse = ICTestAll{TestId}; %
IniTest = ICUse;


IniTest.S(:) = 2.5;
IniTest.C(:) = 8;
IniTest.I(:) = 18;


EquEpoch = 100;
EpocTest = 0;
Contruse = Contr * ones(1600, 1);
Orientationuse = Angle * ones(1600, 1);
responseCtx = GA_h96_context_from_base();

% This is the fixed point response map r(u,p), we only consider u now.
L6FixedLibIndForIteration = [];
responseInfoBaseline = struct();
if L6TonicScaleActive
    [LDEfixedpointBaseline, responseInfoBaseline] = GA_h96_response_map_fixedpoint(...
        Angle, log(lgnSF), log(lgnTF), log(Contr), responseCtx, ...
        'InitialState', IniTest, 'IterationEpoch', EquEpoch);

    LDEUseForFixedL6 = LDEfixedpointBaseline;
    if Isaturation
        LDEUseForFixedL6.E = LDEUseForFixedL6.S * (1-0.3077) + LDEUseForFixedL6.C * 0.3077;
        LDEUseEAdjForFixedL6 = L6Convert(LDEUseForFixedL6.E, EKpUse) ./ LDEUseForFixedL6.E;
        LDEUseForFixedL6.S = LDEUseForFixedL6.S .* LDEUseEAdjForFixedL6;
        LDEUseForFixedL6.C = LDEUseForFixedL6.C .* LDEUseEAdjForFixedL6;
        LDEUseForFixedL6.I = InhMulp(LDEUseForFixedL6.I, IKpUse);
    end
    LDEUseForFixedL6.E = LDEUseForFixedL6.S * (1-0.3077) + LDEUseForFixedL6.C * 0.3077;
    FieldRowForFixedL6 = N_HCOutY * NPixY;
    FieldColForFixedL6 = floor(length(LDEUseForFixedL6.E) / FieldRowForFixedL6);
    L4EfieldForFixedL6 = reshape(LDEUseForFixedL6.E, FieldRowForFixedL6, FieldColForFixedL6);
    L4EfieldPaddedForFixedL6 = padarray(L4EfieldForFixedL6, [1, 1], 'circular');
    CForFixedL6 = conv2(L4EfieldPaddedForFixedL6, L6Kernel, 'same');
    L6FixedUse = L6Convert(CForFixedL6(2:end-1, 2:end-1), L6pars{L6parId});
    L6FixedLibIndForIteration = L6TonicScale * L6FixedUse(:) / 3;
    L6FixedLibIndForIteration(L6FixedLibIndForIteration < 1) = 1;
    L6FixedLibIndForIteration(L6FixedLibIndForIteration > 40) = 40;

    [LDEfixedpoint, responseInfo] = GA_h96_response_map_fixedpoint(...
        Angle, log(lgnSF), log(lgnTF), log(Contr), responseCtx, ...
        'InitialState', IniTest, 'IterationEpoch', EquEpoch, ...
        'FixedL6LibInd', L6FixedLibIndForIteration);
else
    [LDEfixedpoint, responseInfo] = GA_h96_response_map_fixedpoint(...
        Angle, log(lgnSF), log(lgnTF), log(Contr), responseCtx, ...
        'InitialState', IniTest, 'IterationEpoch', EquEpoch);
end


NANFlag(TestId) = responseInfo.NANFlag;
LDEOutAll{TestId} = {LDEfixedpoint};
LDEIL2Diff{TestId} = [];
L2DiffNormNeib{TestId} = [];

LDEOut{TestId} = LDEOutAll{TestId};
LDEEquv = LDEOut{TestId}{end}; % this record the final NESS state
if any(strcmpi(strtrim(getenv('L6NS_SETUP_ONLY')),{'1','true','yes','on'}))
    fprintf('L6NS setup-only mode: fixed point and canonical model variables are ready.\n');
    return
end


if ~NANFlag(TestId)
    for EpcInd = 1:EpocTest+1
        LDEPertRslt = LDEOutAll{TestId}{EpcInd};
        LDEEquv = LDEOut{TestId}{end};
        SDiff = LDEPertRslt.S - LDEEquv.S;
        SDiffHC = reshape(SDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
        SDiff = reshape(SDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
        CDiff = LDEPertRslt.C - LDEEquv.C;
        CDiffHC = reshape(CDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
        CDiff = reshape(CDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
        IDiff = LDEPertRslt.I - LDEEquv.I;
        IDiffHC = reshape(IDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
        IDiff = reshape(IDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);

        L2Diff_OneStepAll(TestId,EpcInd) = ...
            sqrt(sum((SDiff*(1-CplxR)+CDiff*CplxR).^2,'all'));%sqrt(sum([SDiff;CDiff;IDiff].^2));
        DiffVecAll(TestId,EpcInd,:) = [SDiff;CDiff;IDiff];

        aE = 32 / (32 + 8);   % 0.8
        aI = 8  / (32 + 8);   % 0.2
        wC = 0.3077; 
        err_ode_true_all(TestId,EpcInd) = HC_norm_diff(LDEEquv.S, LDEEquv.C, LDEEquv.I, ...
                                LDEPertRslt.S, LDEPertRslt.C, LDEPertRslt.I, ...
                                wC, aE, aI);


    end
end

RunTime = toc;
sprintf('RunTime = %.2f',RunTime)



N_HCOut_full = 4;                   % full network
LDEEquv_rot = HCRot(LDEOutAll{TestId}{end}, RotInd, N_HCOut_full, NPixX, NPixY, MirInd);

S_map_full = reshape(LDEEquv_rot.S, NPixX * N_HCOut_full, NPixY * N_HCOut_full);
C_map_full = reshape(LDEEquv_rot.C, NPixX * N_HCOut_full, NPixY * N_HCOut_full);
I_map_full = reshape(LDEEquv_rot.I, NPixX * N_HCOut_full, NPixY * N_HCOut_full);





    h = 1e-3;
    p = 0.33;
    EpocTest = 1;
    N_HCOut = 4;
    f = [LDEEquv.S; LDEEquv.C; LDEEquv.I];
    PixNumOut = N_HCOut^2 * NPixX * NPixY;
    N         = PixNumOut;

    J = zeros(3*N, 3*N);



    J = compute_J_h96baseline_pref6D_phi(LDEfixedpoint, Contruse,Orientationuse, PixLGNCtgr, L6Kernel, L6pars{L6parId}, ...
    C_SS_meanU, C_CS_meanU, C_IS_mean, ...
    C_SC_meanU, C_CC_meanU, C_IC_mean, ...
    C_SI_mean, C_CI_mean, C_II_mean, ...
    L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
    N_HCOutY, NPixX, NPixY, p, EKpUse, IKpUse,Isaturation,L6EquBlendWeight,L6FixedLibIndForIteration);

    [eig_vec, D] = eig(full(J));


    % Extract eigenvalues as a vector
    eig_J = diag(D);

    [~, idx] = sort(real(eig_J), 'descend');

    eig_J   = eig_J(idx);     % sorted eigenvalues
    eig_vec = eig_vec(:,idx); % sorted eigenvectors

    fprintf('max Re(λ) = %.6f\n', max(real(eig_J)));
    

        % Section 4.3/4.4 task derivatives of the NN drive at fixed state.
    [dPhi, FeatureNames, TaskDerivativeInfo] = GA_h96_pref6D_task_dphi(LDEfixedpoint, Contruse, Orientationuse, ...
        PixLGNCtgr, L6Kernel, L6pars{L6parId}, ...
        C_SS_meanU, C_CS_meanU, C_IS_mean, ...
        C_SC_meanU, C_CC_meanU, C_IC_mean, ...
        C_SI_mean, C_CI_mean, C_II_mean, ...
        L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
        N_HCOutY, NPixX, NPixY, EKpUse, IKpUse, Isaturation);

    % Geometry analysis from proposal sections 4 and 5.
    Section4 = GA_section4_from_jacobian(J, 'dPhi', dPhi, 'FeatureNames', FeatureNames, 'ComputeFullResolvent', true);
    Section4.Transient = GA_discrete_transient_metrics(J, Section4.TaskTangents, FeatureNames, ...
        'MaxStep', 20, 'Threshold', 1.05, 'SampleRank', 16, 'Seed', 1);

    Section5 = GA_section5_from_response(LDEfixedpoint, 'A', J, 'Dur', Section4.TaskTangents, ...
        'CovarianceMode', 'linearized_continuous', ...
        'SolveCovariance', true, 'MaxFullLyapN', 4800);
    Section5.GeometryFeatureNames = FeatureNames;
    if isfield(Section5, 'Geometry') && isstruct(Section5.Geometry) && isfield(Section5.Geometry, 'g') && ~isempty(Section5.Geometry.g)
        Section9 = GA_performance_map(Section4, 'TaskMetric', Section5.Geometry.g);
    else
        Section9 = GA_performance_map(Section4);
    end

    Nx = N_HCOut * NPixX;
    Ny = N_HCOut * NPixY;
    S_mapHCt = reshape(LDEfixedpoint.S, N_HCOut * NPixY, N_HCOut * NPixX);
    C_mapHCt = reshape(LDEfixedpoint.C, N_HCOut * NPixY, N_HCOut * NPixX);
    I_mapHCt = reshape(LDEfixedpoint.I, N_HCOut * NPixY, N_HCOut * NPixX);
    wC = 0.3077;
    E_mapHCt = (1 - wC) * S_mapHCt + wC * C_mapHCt;
    E_map_full = (1 - wC) * S_map_full + wC * C_map_full;

    EigHistEdges = -3:0.05:1;
    EigHistCounts = histcounts(real(eig_J), EigHistEdges);
    TopModeIndex = 1;
    TopModeLambda = eig_J(TopModeIndex);
    TopModeVec = eig_vec(:,TopModeIndex);
    TopModeS = reshape(real(TopModeVec(1:PixNumOut)), Nx, Ny);
    TopModeC = reshape(real(TopModeVec(PixNumOut + (1:PixNumOut))), Nx, Ny);
    TopModeI = reshape(real(TopModeVec(2*PixNumOut + (1:PixNumOut))), Nx, Ny);
    TopModeE = (1 - wC) * TopModeS + wC * TopModeC;

    GeometryMetrics = struct();
    GeometryMetrics.Condition = struct();
    GeometryMetrics.Condition.Contrast = Contr;
    GeometryMetrics.Condition.Angle = Angle;
    GeometryMetrics.Condition.AngleInput = AngleInpt;
    GeometryMetrics.Condition.AngleSource = AngSource;
    GeometryMetrics.Condition.AngleFunctionCategory = AngFuncCtgr;
    GeometryMetrics.Condition.RotInd = RotInd;
    GeometryMetrics.Condition.MirInd = MirInd;
    GeometryMetrics.Condition.lgnSF = lgnSF;
    GeometryMetrics.Condition.lgnTF = lgnTF;
    GeometryMetrics.Condition.xE = xEAll(xEInd);
    GeometryMetrics.Condition.xI = xIAll(xIInd);
    GeometryMetrics.Condition.MonocuFlag = MonocuFlag;
    GeometryMetrics.Parameters = struct();
    GeometryMetrics.Parameters.IKpUse = IKpUse;
    GeometryMetrics.Parameters.EKpUse = EKpUse;
    GeometryMetrics.Parameters.L6pars = L6pars;
    GeometryMetrics.Parameters.L6parId = L6parId;
    GeometryMetrics.Parameters.L6Comment = L6Comment;
    GeometryMetrics.Parameters.L6BaselineChoice = 'c1smooth_LowHighQuadr5_y0_34p0';
    GeometryMetrics.Parameters.L6EquBlendWeight = L6EquBlendWeight;
GeometryMetrics.Parameters.L6DynamicGain = 1 - L6EquBlendWeight;
if numel(L6EquBlendWeight) == 1
    GeometryMetrics.Parameters.L6DerivativeClampMode = 'all_cells_scalar';
else
    GeometryMetrics.Parameters.L6DerivativeClampMode = 'cell_type_specific_S_C_I';
end
GeometryMetrics.Parameters.L6AmplitudeScale = L6AmplitudeScale;
GeometryMetrics.Parameters.L6TonicScale = L6TonicScale;
GeometryMetrics.Parameters.L6TonicScaleActive = L6TonicScaleActive;
GeometryMetrics.Parameters.L6TonicFixedInputUsed = ~isempty(L6FixedLibIndForIteration);
GeometryMetrics.Parameters.L6KernelSigmaScale = L6KernelSigmaScale;
    GeometryMetrics.FixedPoint = GA_fixed_point_summary(LDEfixedpoint);
    GeometryMetrics.FixedPointMaps = struct('S', S_mapHCt, 'C', C_mapHCt, 'I', I_mapHCt, 'E', E_mapHCt);
    GeometryMetrics.RotatedMaps = struct('S', S_map_full, 'C', C_map_full, 'I', I_map_full, 'E', E_map_full);
    GeometryMetrics.ResponseInfo = responseInfo;
    GeometryMetrics.Jacobian = struct();
    GeometryMetrics.Jacobian.MaxRealLambda = max(real(eig_J));
    GeometryMetrics.Jacobian.SpectralRadius = max(abs(eig_J));
    GeometryMetrics.Jacobian.StabilityMargin = 1 - max(real(eig_J));
    GeometryMetrics.Jacobian.NonNormality = Section4.NonNormality;
    GeometryMetrics.Jacobian.EigenValues = eig_J;
    GeometryMetrics.Jacobian.RealPartHistogramEdges = EigHistEdges;
    GeometryMetrics.Jacobian.RealPartHistogramCounts = EigHistCounts;
    GeometryMetrics.TopEigenMode = struct();
    GeometryMetrics.TopEigenMode.Index = TopModeIndex;
    GeometryMetrics.TopEigenMode.Lambda = TopModeLambda;
    GeometryMetrics.TopEigenMode.S = TopModeS;
    GeometryMetrics.TopEigenMode.C = TopModeC;
    GeometryMetrics.TopEigenMode.I = TopModeI;
    GeometryMetrics.TopEigenMode.E = TopModeE;
    GeometryMetrics.TaskDerivative = struct();
    GeometryMetrics.TaskDerivative.dPhi = dPhi;
    GeometryMetrics.TaskDerivative.FeatureNames = FeatureNames;
    GeometryMetrics.TaskDerivative.Info = TaskDerivativeInfo;
    GeometryMetrics.Covariance = struct();
    GeometryMetrics.Covariance.SourceCovarianceModel = Section5.SourceCovarianceModel;
    GeometryMetrics.Covariance.CovarianceModel = Section5.CovarianceModel;
    GeometryMetrics.Covariance.CovarianceEquation = Section5.CovarianceEquation;
    GeometryMetrics.Covariance.CovarianceSolution = Section5.CovarianceSolution;
    GeometryMetrics.Covariance.CovarianceSolved = Section5.CovarianceSolved;
    GeometryMetrics.Covariance.CovarianceConverged = Section5.CovarianceConverged;
    GeometryMetrics.Covariance.NeumannInfo = Section5.NeumannInfo;
    if isfield(Section5, 'Geometry') && isstruct(Section5.Geometry) && isfield(Section5.Geometry, 'g') && ~isempty(Section5.Geometry.g)
        GeometryMetrics.TaskMetric = GA_task_metric_summary(Section5.Geometry.g);
    else
        GeometryMetrics.TaskMetric = [];
    end
    GeometryMetrics.JacobianInputs = struct();
    GeometryMetrics.JacobianInputs.Contruse = Contruse;
    GeometryMetrics.JacobianInputs.Orientationuse = Orientationuse;
    GeometryMetrics.JacobianInputs.PixLGNCtgr = PixLGNCtgr;
    GeometryMetrics.JacobianInputs.L6Kernel = L6Kernel;
    GeometryMetrics.Transient = Section4.Transient;
    GeometryMetrics.FigureMetrics = struct();
    GeometryMetrics.FigureMetrics.EigenSpectrumReal = real(eig_J);
    GeometryMetrics.FigureMetrics.EigenSpectrumImag = imag(eig_J);
    GeometryMetrics.FigureMetrics.EigenRealHistogramEdges = EigHistEdges;
    GeometryMetrics.FigureMetrics.EigenRealHistogramCounts = EigHistCounts;

    GeometryMetrics.Jacobian.EigenVectors = [];
    GeometryMetrics.SavePolicy = 'deduplicated_full_sections_no_standalone_J_no_standalone_eigvec';

    GeometrySaveFile = fullfile(GeometryOutputRoot, ...
        sprintf('geometry_sections4_5_h96baseline_%s_contr%d_angle_%.2f.mat', L6EquComment, Contr, Angle));
    save(GeometrySaveFile, 'Section4', 'Section5', 'Section9', 'GeometryMetrics', ...
        'Contr', 'Angle', '-v7.3');
    fprintf('Saved deduplicated full section 4/5 geometry to %s\n', GeometrySaveFile);


% Eigenvector Heatmap
    
    
    % Nx = N_HCOut * NPixX;   % 4 * 10 = 40
    % Ny = N_HCOut * NPixY;   % 4 * 10 = 40
    % 
    % 
    % S_mapHCt = reshape(LDEfixedpoint.S, N_HCOut * NPixY, N_HCOut * NPixX);
    % C_mapHCt = reshape(LDEfixedpoint.C, N_HCOut * NPixY, N_HCOut * NPixX);
    % I_mapHCt = reshape(LDEfixedpoint.I, N_HCOut * NPixY, N_HCOut * NPixX);
    % wC = 0.3077;
    % E_mapHCt = (1 - wC) * S_mapHCt + wC * C_mapHCt;
    % 
    % 
    % for k = 1:1
    % 
    %     v   = eig_vec(:,k);      % 4800 x 1 eigenvector for mode k
    %     lam = eig_J(k);          % corresponding eigenvalue
    % 
    %     % Split into S, C, I blocks
    %     v_S = v(1:PixNumOut);
    %     v_C = v(PixNumOut + (1:PixNumOut));
    %     v_I = v(2*PixNumOut + (1:PixNumOut));
    % 
    %     % Reshape to 40 x 40 (real part in case of tiny imaginary residue)
    %     V_S_2D = reshape(real(v_S), Nx, Ny);
    %     V_C_2D = reshape(real(v_C), Nx, Ny);
    %     V_I_2D = reshape(real(v_I), Nx, Ny);
    %     V_E_2D = (1 - wC) * V_S_2D + wC * V_C_2D;
    % 
    %     % ---- Plot ----
    %     figure('Visible','off','Color','w','Position',[100 100 2200 900]);
    % 
    % 
    %     subplot(2,4,1);
    %     imagesc(V_S_2D);
    %     axis image;
    %     colorbar;
    %     title(sprintf('Eigenvector %d (S)', k));
    %     hold on;
    %     for h = 1:(N_HCOut-1)
    %         xline(h*NPixX + 0.5, 'k-');
    %         yline(h*NPixY + 0.5, 'k-');
    %     end
    %     hold off;
    % 
    %     subplot(2,4,2);
    %     imagesc(V_C_2D);
    %     axis image;
    %     colorbar;
    %     title(sprintf('Eigenvector %d (C)', k));
    %     hold on;
    %     for h = 1:(N_HCOut-1)
    %         xline(h*NPixX + 0.5, 'k-');
    %         yline(h*NPixY + 0.5, 'k-');
    %     end
    %     hold off;
    % 
    %     subplot(2,4,3);
    %     imagesc(V_I_2D);
    %     axis image;
    %     colorbar;
    %     title(sprintf('Eigenvector %d (I)', k));
    %     hold on;
    %     for h = 1:(N_HCOut-1)
    %         xline(h*NPixX + 0.5, 'k-');
    %         yline(h*NPixY + 0.5, 'k-');
    %     end
    %     hold off;
    % 
    %     subplot(2,4,4);
    %     imagesc(V_E_2D);
    %     axis image;
    %     colorbar;
    %     title(sprintf('Eigenvector %d (E)', k));
    %     hold on;
    %     for h = 1:(N_HCOut-1)
    %         xline(h*NPixX + 0.5, 'k-');
    %         yline(h*NPixY + 0.5, 'k-');
    %     end
    %     hold off;
    % 
    % 
    %     % put fixed point panels below
    % 
    %         % --- S ---
    %     subplot(2,4,5);
    %     imagesc(S_mapHCt);
    %     axis image;
    %     colorbar;
    %     title('S^* (ODE)');
    %     hold on;
    %     % draw HC boundaries
    %     for h = 1:(N_HCOut-1)
    %         xline(h*NPixX + 0.5, 'k-');
    %         yline(h*NPixY + 0.5, 'k-');
    %     end
    %     hold off;
    % 
    %     % --- C ---
    %     subplot(2,4,6);
    %     imagesc(C_mapHCt);
    %     axis image;
    %     colorbar;
    %     title('C^* (ODE)');
    %     hold on;
    %     for h = 1:(N_HCOut-1)
    %         xline(h*NPixX + 0.5, 'k-');
    %         yline(h*NPixY + 0.5, 'k-');
    %     end
    %     hold off;
    % 
    %     % --- I ---
    %     subplot(2,4,7);
    %     imagesc(I_mapHCt);
    %     axis image;
    %     colorbar;
    %     title('I^* (ODE)');
    %     hold on;
    %     for h = 1:(N_HCOut-1)
    %         xline(h*NPixX + 0.5, 'k-');
    %         yline(h*NPixY + 0.5, 'k-');
    %     end
    %     hold off;
    % 
    %     % --- E ---
    %     subplot(2,4,8);
    %     imagesc(E_mapHCt);
    %     axis image;
    %     colorbar;
    %     title('E^* (ODE)');
    %     hold on;
    %     for h = 1:(N_HCOut-1)
    %         xline(h*NPixX + 0.5, 'k-');
    %         yline(h*NPixY + 0.5, 'k-');
    %     end
    %     hold off;
    % 
    % 
    %     sgtitle(sprintf('Eigenvector %d (\lambda = %.4f%+.4fi\n), contr %d, 16 HCs', ...
    %         k, real(lam),imag(lam), Contr));
    % 
    %     filename = sprintf([[repro_paths('project') '/'] ...
    % 'Spectral_Analysis/matlab-inserting_into_CG_model/Complete_Code_for_Paper3/' ...
    % 'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors/' ...
    % 'eig_16HC_5Dmlph96baseline_input_contr%d_angle_%.2f.mat'], ...
    % Contr, Angle);
    %     [saveFolder, ~, ~] = fileparts(filename);
    % 
    %     fig_name_k = sprintf('EigMode%d_16HC_heatmap_5Dmlph96baseline_contr%d_angle%.2f.png', ...
    %         k, Contr, Angle);
    %     fig_path_png = fullfile(saveFolder, fig_name_k);
    %     fig_path_fig = strrep(fig_path_png, '.png', '.fig');
    %     saveas(gcf, fig_path_png);
    %     savefig(gcf, fig_path_fig);
    % 
    % end


toc


% filename = sprintf([[repro_paths('project') '/'] ...
%     'Spectral_Analysis/matlab-inserting_into_CG_model/Complete_Code_for_Paper3/' ...
%     'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors/' ...
%     'eig_16HC_5Dmlph96baseline_input_contr%d_angle_%.2f.mat'], ...
%     Contr, Angle);
% 
% 
%     save(filename, 'eig_J', 'eig_vec');

 

% Full Eigenspectrum


% % Extract folder path from the .mat filename
% [saveFolder, ~, ~] = fileparts(filename);
% 
% % Create figure
% hFig = figure('Visible','off');
% plot(real(eig_J), imag(eig_J), 'o', ...
%      'MarkerSize', 6, ...
%      'LineWidth', 1.2);
% hold on;
% 
% % Reference axes
% yline(0, 'k--');
% xline(0, 'k--');
% 
% axis equal;
% grid on;
% xlim([-3 1]);
% 
% xlabel('Re(\lambda)');
% ylabel('Im(\lambda)');
% title(sprintf('Eigenvalues (Angle = %.2f°, Contr = %d)', Angle, Contr));
% 
% 
% 
% figFilename = fullfile(saveFolder, ...
%     sprintf('eigenvalue_spectrum_5Dmlph96baseline_angle_%.2f_contr%d.png', Angle, Contr));
% 
% saveas(hFig, figFilename);
% 
% % (Optional) also save MATLAB figure
% savefig(hFig, strrep(figFilename, '.png', '.fig'));
% 
% % close(hFig);

% Eigenvalue real part histogram

% hFig2 = figure('Visible','off');
% histogram(real(eig_J), -3:0.05:1);   % 40 bins (adjust as needed)
% xlabel('Re(\lambda)');
% ylabel('Count');
% title(sprintf('Distribution of Re(\\lambda), Angle = %.2f°, Contr = %d', Angle, Contr));
% grid on;
% 
% % Save histogram
% histFilename = fullfile(saveFolder, ...
%     sprintf('eigenvalue_realpart_hist_5Dmlph96baseline_angle_%.2f_contr%d.png', Angle, Contr));
% saveas(hFig2, histFilename);
% savefig(hFig2, strrep(histFilename, '.png', '.fig'));

    
end
end
fprintf('Completed Torch array task %d: angle %.2f, contrast %d\n', TaskId, AngleCond, ContrCond);
