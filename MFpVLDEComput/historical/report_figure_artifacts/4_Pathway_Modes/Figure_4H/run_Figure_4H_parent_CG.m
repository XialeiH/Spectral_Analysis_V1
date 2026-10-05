%% Paper3: 1st test
% How will a 4*8 binocular patch work?
%% Setting up WS

tic

CurrentFolder = getenv('FIG1F1_REPO_ROOT');
if isempty(CurrentFolder)
    CurrentFolder = pwd;
end
helperOverride = getenv('FIG1F1_HELPER_ROOT');
if ~isempty(helperOverride)
    addpath(helperOverride, '-begin');
end
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
SaveToFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper3PlotingData/']; % V1D2
DataFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/25Function_Binocular/']; % V1D2
DataFolder1 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/16Function_Scheme/LARGE16/']; % V1D2
DataFolder2 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/25Function_binocular_realLGN/']; % V1D2
DataFolder3 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper2PlotingData/Typical_trajs/']; % V1D2
DataFolder4 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper3ICTestData/'];
addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'])

addpath(DataFolder4)
addpath(DataFolder3)
addpath(SaveToFolder)
addpath(DataFolder2)
addpath(DataFolder1)
addpath(DataFolder)
addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'])
if ~isempty(helperOverride)
    addpath(helperOverride, '-begin');
end
% B1. Load LDE computations and do Fr distributions
% Should rerun Sect 1 after in case loaded file contain any paths.

DataPt = 'V4D2'; %'V4D1'; 'V5D1'
load(sprintf("AllMFPixPara_Paper2TuneFig1%s.mat",DataPt))

FigurePaperPath = [CurrentFolder '/Figures/Demo062324/']; % V1D2
MonocuFlag = false;
close all
%% 2. Test for different x for E and I; yE and yI kept as 1 -- All killed are folded back
% first get basic function information
% 
% Too time comsuming. --->> using HPC

close all
xEAll = 0:0.2:1; xIAll = 0:0.2:1;
lgnSF = 2.5; lgnTF = 10;

ExpTex = 'CtrlL4';
Comment = ['L6Real' ExpTex];

L6pars = {};
L6ParamRawUse = {2.5,3,62.5,86,80.4,106,90,142.2,102,166.7,108,210,118, ...
    {[0 15 29],[0,34.0],[28 32 50],[69 85]}};
L6C1Grid = 0:0.25:120;
L6pars{5} = {2.5,3,{'c1smooth',L6ParamRawUse,L6C1Grid}};

% figure(101)
% xx = 0:0.1:100;
% hold on
% for L6Id = 1:3
%     yy = L6Convert(xx,L6pars{L6Id});
%     plot(xx,yy,'DisplayName',sprintf('L6Convert-%d',L6Id))
% end
% grid on
% xlim([0 90])
% ylim([0 120])
% xlabel('L4E')
% ylabel('frL6')
% legend('Location','SouthEast')

IKp7.Thrsld1 = 54;
IKp7.Thrsld2 = 63;
IKp7.Highist = 103;
IKp7.Slope = 0.941;
IKp7.down1 = 0.3;
Isaturation = true;
IKp7.Mode = 'multisigmoid';
IKp7.IntaH = 4.5;
IKp7.IntaL = -2.5;
IKp7.IntbL = 4;
IKp7.IntbH = -0.3;
IKp7.SmoothJoinHalfWidth = [1 0 1];
IKp7.SmoothT2QuinticWidth = 3;

IKpUse = IKp7;
%%
close all
% figure(102)
EKpUse = {1,0,[50; 70; 100],[50; 60.5; 69],'quadratic'}; 
% % we are using the L4eL6 conversion curve formula for L4E depression
% xx = 0:0.1:220;
% hold on
%     yy = L6Convert(xx,EKpUse);
%     plot(xx,yy,'DisplayName',sprintf('L4E-depression'))
%     xlabel('L4E')
% ylabel('L4E-dep')
% grid on
%%
close all
%Comment = '';
xEInd = 6;
xIInd = 6;
fprintf('xE=%.1f, xI=%.1f\n',xEAll(xEInd),xIAll(xIInd))
L6parId = 5;
fprintf('L6 ParameterID %d\n',L6parId)
if ischar(L6pars{L6parId}{end}) && strcmp(L6pars{L6parId}{end}, 'quadratic')
    L6Comment = sprintf('%s',L6pars{L6parId}{end})
else
    L6Comment = sprintf('linear')
end
tic
%AngleList = {'0.0'}; % only vertical for now
AngleList = {'0.0','7.5','15.0','22.5'};
DomList =  {'Small','Large','LARGER'};
%% 3. Set up initial conditions

% there are 3200 pixels in total
N_HCOutX = 4; N_HCOutY = 4; PixNumOut = N_HCOutX*N_HCOutY*NPixX*NPixY;
N_HCinX = 4; N_HCinY = 4;
NPixX = 10; NPixY = 10;
%% 
% Parameters of L4 kernels

% here I test multiple angles, 
% For the Sigmoid version, only test 0 degree 
AngleTestAll = 0:7.5:22.5; % 0:7.5:22.5; % has to be multiples of 7.5
TestNum = length(AngleTestAll); ICTestAll = cell(TestNum,1);
TestExSeq = xEAll(xEInd)*ones(size(0:7.5:180)); TestEySeq = 1*ones(size(TestExSeq));
TestIxSeq = xIAll(xIInd)*ones(size(0:7.5:180)); TestIySeq = 1*ones(size(TestExSeq));
%% 
% %%%60 deg: by rotating 15 deg pattern
% 
% Start from 90, then go to 22.5

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
%% 
% load equiv

%PlotUse = load([DataFolder2 sprintf('LDETraces_bg-ang%s.mat',AngPrint)],'LDEEquv');
%LDEEquv = PlotUse.LDEEquv;
%LDEEquv = HCRot(LDEEquv,RotInd,N_HCOut,NPixX,NPixY,MirInd);
%% 
% Hyperparameters of iterations

p = 0.33;% 1-p for the original input
EpocTest = 100;
ExportFlag = 'xn';

L2Diff_OneStepAll = zeros(length(ICTestAll),EpocTest+1);
DiffVecAll = zeros(length(ICTestAll),EpocTest+1,3*NPixY*NPixX);
%% 
% Since we are doing different network architectures now - need to separate 
% inputs to the S C I functions!!!

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
%% Compute NESSs for each angle

% here I only test 0 degree
TestId = 1;
AngleInpt = AngleTestAll(TestId);

RotInd = floor(AngleInpt/45); % 0-3
MirInd = mod(floor(AngleInpt/22.5),2); % 0: no mirror; 1: mirror
switch MirInd % find the corresponding source
    case 0
        AngSource = mod(AngleInpt,45);
    case 1
        AngSource = mod(-AngleInpt,45);
end

% Use the source, decide which response function for angle should I use
AngFuncCtgrCdid = 0:7.5:22.5;
[~,AngFuncCtgr] = min(abs(AngFuncCtgrCdid - AngSource));

AngPrint = AngleList{AngFuncCtgr};
Angle = str2num(AngPrint); % can select from 0 7.5 15 22.5
%% 
% load all precomputed domains

L4EmeshXAll = cell(size(DomList));
L4ImeshYAll = cell(size(DomList));
LDEFrfuncAll = cell(size(DomList));
for DomInd = 1:length(DomList)
    FuncTemp = load(sprintf('Func200%s-%s-Ang%s-SF%.1fTF%d-%s.mat',...
        DataPt,Comment,AngPrint,lgnSF,lgnTF,DomList{DomInd}),...
        'LDEFrfunc','L4EmeshX','L4ImeshY');
    L4EmeshXAll{DomInd} = FuncTemp.L4EmeshX;
    L4ImeshYAll{DomInd} = FuncTemp.L4ImeshY;
    LDEFrfuncAll{DomInd} = FuncTemp.LDEFrfunc;
end
%% 4. Iteration 
% New: Need to extend all these to arbitrary HC numbers!

%    N_HCOut = 4; PixNumOut = N_HCOut^2*NPixX*NPixY;
LGNlist = 5; L6list = 33; % Why do we have L6list here?
OD_SMapModi = OD_SMap;
%% 
% Needs to replace the midel row as "5" -- Background

%OD_SMapModi = OD_SMap(1:2*n_S_HC,1:2*n_S_HC);
if MonocuFlag
    OD_SMapModi(n_S_HC+1:2*n_S_HC,1:end) = 5;
    SaveStr = ['Mono' Comment];
else
    SaveStr = ['Bino' Comment];
end
%% 
% L6 kernel

L6ConvSig = 0.75;
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
%% 
% LGN and L6 are blurred 

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
%% 
% completely BG or not?

BGFlag = false;
if BGFlag
    PixL6Ctgr(:,end) = 1; PixL6Ctgr(:,1:end-1) = 0;
    PixLGNCtgr(:,end) = 1; PixLGNCtgr(:,1:end-1) = 0;
end
%% 
% However, the Ctgr for pixels needs to be rotated/flipped to represent the 
% actual gratings
% 
% *Instead of rotating existing ctgr fields, I should simply permute 1234*  

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
%% 
% Get pixelwise connectivities
% 
% NOTE: NEED to modify AveSpatKer_Rec to implement the ocular modulation of 
% connectivity kernels    

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
%% 
% Firing rate suppresion

%     EKp.Thrsld = 50; EKp.Highist = [200,150]; EKp.HardBound = 200; EKp.Slope = 0;
%
%     IKp.Thrsld = 70; IKp.Highist = [100,95];  IKp.HardBound = 120; IKp.Slope = 1; % 1; 0.9(95); 0.8(9)
%% 
% set up IC

tic
ICUse = ICTestAll{TestId}; %
IniTest = ICUse;
Figure4HParentCG = strcmp(getenv('FIGURE4H_PARENT_CG'), '1');
if Figure4HParentCG
    baselineData = load(getenv('FIGURE4H_CG_BASELINE_FILE'), 'LDEEquv');
    perturbationData = load(getenv('FIGURE4H_PERTURBATION_FILE'), ...
        'perturbationDirection');
    parentCGFixed = baselineData.LDEEquv;
    direction = double(perturbationData.perturbationDirection(:));
    populationSize = numel(direction)/3;
    directionE = 0.6923*direction(1:populationSize) + ...
        0.3077*direction(populationSize+(1:populationSize));
    directionI = direction(2*populationSize+(1:populationSize));
    directionHC = sqrt(0.8*mean(directionE.^2) + 0.2*mean(directionI.^2));
    direction = direction/max(directionHC,eps); % matched HC=1 perturbation
    IniTest = struct( ...
        'S',parentCGFixed.S + direction(1:populationSize), ...
        'C',parentCGFixed.C + direction(populationSize+(1:populationSize)), ...
        'I',parentCGFixed.I + direction(2*populationSize+(1:populationSize)));
    EpocTest = 180;
end
%IniTest.S = symmHCs(ICUse.S,N_HCOut,NPixX,NPixY);
%IniTest.C = symmHCs(ICUse.C,N_HCOut,NPixX,NPixY);
%IniTest.I = symmHCs(ICUse.I,N_HCOut,NPixX,NPixY);

%% 
% start iteration

% Use the actual Paper 3 response tables. The similarly named
% _Sigmoid_moe_v4 function is a surrogate and does not read LDEFrfuncAll.
[LDEOutAll{TestId},LDEIL2Diff{TestId},L2DiffNormNeib{TestId},~,~,FuncUse,...
    NANFlag(TestId)] = ...
    LDEIteration_135FuncMain_CombDom_RealLGNL6(...
    PixLGNCtgr,L6Kernel,IniTest,p,L6pars{L6parId}, EpocTest,...
    C_SS_meanU,C_CS_meanU,C_IS_mean,...
    C_SC_meanU,C_CC_meanU,C_IC_mean,...
    C_SI_mean, C_CI_mean, C_II_mean,...
    L4SEp, L4SIp, ...
    L4CEp, L4CIp, ...
    L4IEp, L4IIp, ...
    L4EmeshXAll,L4ImeshYAll,LDEFrfuncAll,CtgrOrderReadout,...
    N_HCOutY,NPixX,NPixY,Isaturation,'xn',EKpUse,IKpUse)  ;

LDEEquv = LDEOutAll{TestId}{end};
if NANFlag(TestId) || any(~isfinite([LDEEquv.S; LDEEquv.C; LDEEquv.I]))
    error('Figure1F1:InvalidLibraryEquilibrium', ...
        'The genuine library iteration did not produce a finite equilibrium.');
end
fprintf('Library equilibrium ranges: S=[%.6g %.6g], C=[%.6g %.6g], I=[%.6g %.6g]\n', ...
    min(LDEEquv.S), max(LDEEquv.S), min(LDEEquv.C), max(LDEEquv.C), ...
    min(LDEEquv.I), max(LDEEquv.I));

if Figure4HParentCG
    trajectoryCells = LDEOutAll{TestId};
    parentCGStates = zeros(3*populationSize,numel(trajectoryCells));
    for stateIndex = 1:numel(trajectoryCells)
        state = trajectoryCells{stateIndex};
        parentCGStates(:,stateIndex) = [state.S; state.C; state.I];
    end
    parentCGFixedVector = [parentCGFixed.S; parentCGFixed.C; parentCGFixed.I];
    % Saved as iteration indices. The analysis converts each epoch to
    % p*tau milliseconds using the calibrated ODE time constant.
    parentCGTimesMs = 0:(numel(trajectoryCells)-1);
    outputFile = getenv('FIGURE4H_CG_OUTPUT_FILE');
    save(outputFile,'parentCGStates','parentCGTimesMs', ...
        'parentCGFixedVector','direction','IKpUse','L6pars','L6parId','-v7.3');
    fprintf('Saved Figure 4H parent-CG trajectory to %s.\n',outputFile);
    return
end

smokeOnly = strcmp(getenv('FIG1F1_SMOKE_ONLY'), '1');
if smokeOnly
    outputFile = getenv('FIG1F1_OUTPUT_FILE');
    save(outputFile, 'LDEEquv', 'FuncUse', 'IKpUse', 'L6pars', 'L6parId', '-v7.3');
    fprintf('Saved genuine-library equilibrium smoke result to %s\n', outputFile);
    return
end

[J, jacobianDiagnostics] = compute_J_library135_phi(LDEEquv, PixLGNCtgr, ...
    L6Kernel, L6pars{L6parId}, ...
    C_SS_meanU, C_CS_meanU, C_IS_mean, ...
    C_SC_meanU, C_CC_meanU, C_IC_mean, ...
    C_SI_mean, C_CI_mean, C_II_mean, ...
    L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
    L4EmeshXAll, L4ImeshYAll, LDEFrfuncAll, CtgrOrderReadout, ...
    N_HCOutY, NPixX, NPixY, EKpUse, IKpUse, Isaturation);

clusterCompare = strcmp(getenv('FIG1F1_CLUSTER_COMPARE'), '1');
if clusterCompare
    clusterModeCount = 10;
else
    clusterModeCount = 4;
end
clusterOnly = strcmp(getenv('FIG1F1_CLUSTER4_ONLY'), '1') || clusterCompare;
if clusterOnly
    eigOptions = struct('tol', 1e-10, 'maxit', 5000, ...
        'p', 80, 'disp', 0, 'isreal', true);
    [eigVectors, eigDiagonal, eigFlag] = eigs( ...
        sparse(J), clusterModeCount, 'largestreal', eigOptions);
    if eigFlag ~= 0
        error('Figure1F1:Eigs', 'eigs returned flag %d.', eigFlag);
    end
else
    [eigVectors, eigDiagonal] = eig(full(J));
end
eigValues = diag(eigDiagonal);
[~, eigOrder] = sortrows([real(eigValues), imag(eigValues)], [-1 -2]);
eigValues = eigValues(eigOrder);
eigVectors = eigVectors(:, eigOrder);

populationSize = numel(LDEEquv.S);
mapSide = N_HCOutY * NPixY;
wC = 0.3077;
wS = 1 - wC;
fixedPointMaps = struct( ...
    'S', reshape(real(LDEEquv.S), mapSide, mapSide), ...
    'C', reshape(real(LDEEquv.C), mapSide, mapSide), ...
    'I', reshape(real(LDEEquv.I), mapSide, mapSide), ...
    'E', reshape(real(wS * LDEEquv.S + wC * LDEEquv.C), mapSide, mapSide));

clusterVectors = eigVectors(:, 1:clusterModeCount);
for modeIndex = 1:clusterModeCount
    clusterVectors(:, modeIndex) = clusterVectors(:, modeIndex) / ...
        max(norm(clusterVectors(:, modeIndex)), eps);
end
[clusterBasis, ~] = qr(clusterVectors, 0);
clusterE = wS * clusterBasis(1:populationSize, :) + ...
    wC * clusterBasis(populationSize + (1:populationSize), :);
topEigenclusterE = reshape(sqrt(sum(abs(clusterE).^2, 2)), mapSide, mapSide);

if clusterOnly
    outputFile = getenv('FIG1F1_OUTPUT_FILE');
    if clusterCompare
        clusterCounts = [6 8 10];
        topEigenclusterEAll = cell(size(clusterCounts));
        for clusterIndex = 1:numel(clusterCounts)
            count = clusterCounts(clusterIndex);
            vectors = eigVectors(:, 1:count);
            for modeIndex = 1:count
                vectors(:, modeIndex) = vectors(:, modeIndex) / ...
                    max(norm(vectors(:, modeIndex)), eps);
            end
            [basis, ~] = qr(vectors, 0);
            eBasis = wS * basis(1:populationSize, :) + ...
                wC * basis(populationSize + (1:populationSize), :);
            topEigenclusterEAll{clusterIndex} = reshape( ...
                sqrt(sum(abs(eBasis).^2, 2)), mapSide, mapSide);
        end
        save(outputFile, 'topEigenclusterEAll', 'clusterCounts', ...
            'eigValues', '-v7');
        fprintf('Saved 6/8/10-mode genuine-library eigenclusters to %s\n', outputFile);
    else
        save(outputFile, 'topEigenclusterE', 'clusterModeCount', 'eigValues', '-v7');
        fprintf('Saved four-mode genuine-library eigencluster to %s\n', outputFile);
    end
    return
end

svdOptions = struct('tol', 1e-9, 'maxit', 3500, ...
    'p', min(size(J, 1), 80), 'disp', 0);
[singularOutput, singularValue, ~, singularFlag] = ...
    svds(sparse(J), 1, 'largest', svdOptions);
if singularFlag ~= 0
    error('Figure1F1:Svds', 'svds returned flag %d.', singularFlag);
end
[~, singularPivot] = max(abs(singularOutput));
singularOutput = singularOutput * exp(-1i * angle(singularOutput(singularPivot)));
if real(singularOutput(singularPivot)) < 0
    singularOutput = -singularOutput;
end
singularE = wS * singularOutput(1:populationSize) + ...
    wC * singularOutput(populationSize + (1:populationSize));
topSingularModeE = reshape(real(singularE), mapSide, mapSide);

outputFile = getenv('FIG1F1_OUTPUT_FILE');
if isempty(outputFile)
    outputFile = fullfile(pwd, 'figure_1f1_current_3d_library.mat');
end
save(outputFile, 'LDEEquv', 'fixedPointMaps', 'eigValues', ...
    'topEigenclusterE', 'topSingularModeE', 'clusterModeCount', ...
    'singularValue', 'IKpUse', 'L6pars', 'L6parId', ...
    'jacobianDiagnostics', 'FuncUse', '-v7.3');
fprintf('Saved current-parameter 3D library result to %s\n', outputFile);
fprintf('max Re(lambda) = %.12g; sigma_1 = %.12g\n', ...
    max(real(eigValues)), singularValue);
return
%%

% [LDEOut{TestId},LDEIL2Diff{TestId},L2DiffNormNeib{TestId},~,~,FuncUse,...
%     NANFlag(TestId)] = ...
%     LDEIteration_135FuncMain_CombDom_RealLGNL6(...
%     PixLGNCtgr,L6Kernel,IniTest,p,L6pars{L6parId}, EpocTest,...
%     C_SS_meanU,C_CS_meanU,C_IS_mean,...
%     C_SC_meanU,C_CC_meanU,C_IC_mean,...
%     C_SI_mean, C_CI_mean, C_II_mean,...
%     L4SEp, L4SIp, ...
%     L4CEp, L4CIp, ...
%     L4IEp, L4IIp, ...
%     L4EmeshXAll,L4ImeshYAll,LDEFrfuncAll,...
%     N_HCOutY,NPixX,NPixY,Isaturation,'xn',EKpUse,IKpUse)  ;
% 
%     LDEEquv = LDEOut{TestId}{end}; % this record the final NESS state
%     % LDEEquv = ICUse;  % this is to compare it with the precomputed NESS state, that is ZC's result.
% 
% RunTime = toc;
% sprintf('RunTime = %.2f',RunTime)
% 
% % this is for computing the distance from each epoch to the Final NESS state
% if ~NANFlag(TestId)
%     for EpcInd = 1:EpocTest+1
%         LDEPertRslt = LDEOutAll{TestId}{EpcInd};
%         SDiff = LDEPertRslt.S - LDEEquv.S;
%         SDiffHC = reshape(SDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
%         SDiff = reshape(SDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
%         CDiff = LDEPertRslt.C - LDEEquv.C;
%         CDiffHC = reshape(CDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
%         CDiff = reshape(CDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
%         IDiff = LDEPertRslt.I - LDEEquv.I;
%         IDiffHC = reshape(IDiff,N_HCOutY*NPixY,N_HCOutX*NPixX);
%         IDiff = reshape(IDiffHC(1:NPixY,1:NPixX),NPixY*NPixX,1);
% 
%         L2Diff_OneStepAll(TestId,EpcInd) = ...
%             sqrt(sum((SDiff*(1-CplxR)+CDiff*CplxR).^2,'all'));%sqrt(sum([SDiff;CDiff;IDiff].^2));
%         DiffVecAll(TestId,EpcInd,:) = [SDiff;CDiff;IDiff];
% 
%         aE = 32 / (32 + 8);   % 0.8
%         aI = 8  / (32 + 8);   % 0.2
%         wC = 0.3077; 
%         err_ode_true_all(TestId,EpcInd) = HC_norm_diff(LDEEquv.S, LDEEquv.C, LDEEquv.I, ...
%                                 LDEPertRslt.S, LDEPertRslt.C, LDEPertRslt.I, ...
%                                 wC, aE, aI);
% 
% 
%     end
% end
%%
N_HCOut_full = 4;                   % full network
LDEEquv_rot = HCRot(LDEEquv, RotInd, N_HCOut_full, NPixX, NPixY, MirInd);


S_map_full = reshape(LDEEquv_rot.S, NPixX * N_HCOut_full, NPixY * N_HCOut_full);
C_map_full = reshape(LDEEquv_rot.C, NPixX * N_HCOut_full, NPixY * N_HCOut_full);
I_map_full = reshape(LDEEquv_rot.I, NPixX * N_HCOut_full, NPixY * N_HCOut_full);

% Plot
figure;

subplot(1,3,1);
imagesc(S_map_full); colorbar; axis equal tight;
title('4×4 HCs (S)');

subplot(1,3,2);
imagesc(C_map_full); colorbar; axis equal tight;
title('4×4 HCs (C)');

subplot(1,3,3);
imagesc(I_map_full); colorbar; axis equal tight;
title('4×4 HCs (I)');
%%


    h = 1e-3;
    p = 0.33;
    EpocTest = 1;
    N_HCOut = 4;
    f = [LDEEquv.S; LDEEquv.C; LDEEquv.I];
    PixNumOut = N_HCOut^2 * NPixX * NPixY;
    N         = PixNumOut;

    J = zeros(3*N, 3*N);
%%
J = compute_J_moe_v4_phi(LDEOutAll{TestId}{end}, PixLGNCtgr, L6Kernel, L6pars{L6parId}, ...
    C_SS_meanU, C_CS_meanU, C_IS_mean, ...
    C_SC_meanU, C_CC_meanU, C_IC_mean, ...
    C_SI_mean, C_CI_mean, C_II_mean, ...
    L4SEp, L4SIp, L4CEp, L4CIp, L4IEp, L4IIp, ...
    N_HCOutY, NPixX, NPixY, p, EKpUse, IKpUse,Isaturation);


     [eig_vec, D] = eig(full(J));

    % Extract eigenvalues as a vector
    eig_J = diag(D);

    [~, idx] = sort(real(eig_J), 'descend');

    eig_J   = eig_J(idx);     % sorted eigenvalues
    eig_vec = eig_vec(:,idx); % sorted eigenvectors

    fprintf('max Re(λ) = %.6f\n', max(real(eig_J)));

    
%%
    Nx = N_HCOut * NPixX;   % 4 * 10 = 40
    Ny = N_HCOut * NPixY;   % 4 * 10 = 40


    S_mapHCt = reshape(LDEEquv.S, N_HCOut * NPixY, N_HCOut * NPixX);
    C_mapHCt = reshape(LDEEquv.C, N_HCOut * NPixY, N_HCOut * NPixX);
    I_mapHCt = reshape(LDEEquv.I, N_HCOut * NPixY, N_HCOut * NPixX);


    for k = 1:10

        v   = eig_vec(:,k);      % 4800 × 1 eigenvector for mode k
        lam = eig_J(k);          % corresponding eigenvalue

        % Split into S, C, I blocks
        v_S = v(1:PixNumOut);
        v_C = v(PixNumOut + (1:PixNumOut));
        v_I = v(2*PixNumOut + (1:PixNumOut));

        % Reshape to 40 × 40 (real part in case of tiny imaginary residue)
        V_S_2D = reshape(real(v_S), Nx, Ny);
        V_C_2D = reshape(real(v_C), Nx, Ny);
        V_I_2D = reshape(real(v_I), Nx, Ny);

        % ---- Plot ----
        figure('Color','w','Position',[100 100 1800 900]);


        subplot(2,3,1);
        imagesc(V_S_2D);
        axis image;
        colorbar;
        title(sprintf('Eigenvector %d (S)', k));
        hold on;
        for h = 1:(N_HCOut-1)
            xline(h*NPixX + 0.5, 'k-');
            yline(h*NPixY + 0.5, 'k-');
        end
        hold off;

        subplot(2,3,2);
        imagesc(V_C_2D);
        axis image;
        colorbar;
        title(sprintf('Eigenvector %d (C)', k));
        hold on;
        for h = 1:(N_HCOut-1)
            xline(h*NPixX + 0.5, 'k-');
            yline(h*NPixY + 0.5, 'k-');
        end
        hold off;

        subplot(2,3,3);
        imagesc(V_I_2D);
        axis image;
        colorbar;
        title(sprintf('Eigenvector %d (I)', k));
        hold on;
        for h = 1:(N_HCOut-1)
            xline(h*NPixX + 0.5, 'k-');
            yline(h*NPixY + 0.5, 'k-');
        end
        hold off;


        % put fixed point penals below

            % --- S ---
        subplot(2,3,4);
        imagesc(S_mapHCt);
        axis image;
        colorbar;
        title('S^* (ODE)');
        hold on;
        % draw HC boundaries
        for h = 1:(N_HCOut-1)
            xline(h*NPixX + 0.5, 'k-');
            yline(h*NPixY + 0.5, 'k-');
        end
        hold off;

        % --- C ---
        subplot(2,3,5);
        imagesc(C_mapHCt);
        axis image;
        colorbar;
        title('C^* (ODE)');
        hold on;
        for h = 1:(N_HCOut-1)
            xline(h*NPixX + 0.5, 'k-');
            yline(h*NPixY + 0.5, 'k-');
        end
        hold off;

        % --- I ---
        subplot(2,3,6);
        imagesc(I_mapHCt);
        axis image;
        colorbar;
        title('I^* (ODE)');
        hold on;
        for h = 1:(N_HCOut-1)
            xline(h*NPixX + 0.5, 'k-');
            yline(h*NPixY + 0.5, 'k-');
        end
        hold off;


        sgtitle(sprintf('Eigenvector %d (\\lambda = %.4f%+.4fi\n), 16 HCs', ...
            k, real(lam),imag(lam)));

        filename = sprintf(['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis/matlab-inserting_into_CG_model/Complete_Code_for_Paper3/' ...
    'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors/3DSigmoid_numerical_results/EigMode_16HC_heatmap_3DSigmoid_numerical_angle0.0' ...
    'eig_16HC_3DLibrary_input_angle_%.1f.mat'], ...
    Angle);
        [saveFolder, ~, ~] = fileparts(filename);
     
        fig_name_k = sprintf('EigMode%d_16HC_heatmap_3DLibrary_angle%.1f.png', ...
            k, Angle);
        fig_path_png = fullfile(saveFolder, fig_name_k);
        fig_path_fig = strrep(fig_path_png, '.png', '.fig');
        saveas(gcf, fig_path_png);
        savefig(gcf, fig_path_fig);
    end
toc
%%
      %       fig_name_k = sprintf('EigMode%d_16HC_heatmap_3DLibraey_input%.1f_tau%.4f_EE%.4f_II%.4f_EI%.4f_IE%.4f.png', ...
      %       k, Angle, tau, EEconnectivity_scale, IIconnectivity_scale, ...
      %       EIconnectivity_scale, IEconnectivity_scale);
      % 
      %   print(gcf, fullfile(fig_folder, fig_name_k), ...
      % '-dpng', '-r600');  
      %   fprintf('Saved eigenvector heatmap (mode %d): %s\n', ...
      %       k, fullfile(fig_folder, fig_name_k));
%%
filename = sprintf(['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis/matlab-inserting_into_CG_model/Complete_Code_for_Paper3/' ...
    'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors/' ...
    'eig_16HC_3DSigmoid_numerical_input_%.1f.mat'], ...
    Angle);
     

    save(filename, 'eig_J', 'eig_vec');

 
%%
Angle = 22.5;
    filename = sprintf(['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis/matlab-inserting_into_CG_model/Complete_Code_for_Paper3/' ...
    'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors/' ...
    'eig_16HC_3DLibrary_input_%.1f.mat'], ...
    Angle);
% Extract folder path from the .mat filename
[saveFolder, ~, ~] = fileparts(filename);

% Create figure
hFig = figure;
plot(real(eig_J), imag(eig_J), 'o', ...
     'MarkerSize', 6, ...
     'LineWidth', 1.2);
hold on;

% Reference axes
yline(0, 'k--');
xline(0, 'k--');

axis equal;
grid on;
xlim([-3 1]);
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
title(sprintf('Eigenvalues (Angle = %.1f°)', Angle));



figFilename = fullfile(saveFolder, ...
    sprintf('eigenvalue_spectrum_angle_%.1f.png', Angle));

saveas(hFig, figFilename);

% (Optional) also save MATLAB figure
savefig(hFig, strrep(figFilename, '.png', '.fig'));

% close(hFig);
%%
filename = sprintf(['/Users/xialeihuang/Desktop/Neuroscience_Project/' ...
    'Spectral_Analysis/matlab-inserting_into_CG_model/Complete_Code_for_Paper3/' ...
    'NYU-Vision-2Drive-main/Figures/spectral_analysis_eigenvalue_eigenvectors/' ...
    'eig_16HC_3DSigmoid_numerical_input_%.1f.mat'], ...
    Angle);
[saveFolder, ~, ~] = fileparts(filename);

hFig2 = figure;
histogram(real(eig_J), -3:0.05:1);   % 40 bins (adjust as needed)
xlabel('Re(\lambda)');
ylabel('Count');
title(sprintf('Distribution of Re(\\lambda), Angle = %.1f°', Angle));
grid on;

% Save histogram
histFilename = fullfile(saveFolder, ...
    sprintf('eigenvalue_realpart_hist_angle_%.1f.png', Angle));
saveas(hFig2, histFilename);
savefig(hFig2, strrep(histFilename, '.png', '.fig'));
