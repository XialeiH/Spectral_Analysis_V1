%% Paper 2: Fig 9
%% 1. Setting up network
%% Setting up WS

CurrentFolder = pwd
addpath(CurrentFolder)
addpath([CurrentFolder '/Utils'])
addpath([CurrentFolder '/Data'])
DataFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper2PlotingData/']; % V1D2
DataFolder1 = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/GlobConv/']; % V1D2

addpath(DataFolder1)
addpath(DataFolder)
addpath([CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/'])

FigurePaperPath = [CurrentFolder '/Paper2Figs/']; % V1D2
SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper2NWSimulationData/'];
%% 
% Then use the new function to generate new network stucture


%% If now, we start computation below...
% Frst setup network
N_HC = 3;
% Number of E and I neurons
n_S_HC = 45; n_C_HC = 30; n_I_HC = 31; % per side of HC
N_S = n_S_HC^2 * N_HC^2; % neuron numbers In all
N_C = n_C_HC^2 * N_HC^2; % neuron numbers In all
N_I = n_I_HC^2 * N_HC^2;
N_E = N_S+N_C; CplxR = N_C/N_E;
% Grid sizes of E and I neurons;
Size_HC = 0.500; % in mm;
Size_S = Size_HC/n_S_HC; Size_C = Size_HC/n_C_HC; Size_I = Size_HC/n_I_HC;
% Projection: SD of distances
SD_E = 0.2/sqrt(2); SD_I = 0.125/sqrt(2);
Dist_LB = 0.36; % ignore the connection probability of dist>0.3mm
% Peak probability of projection
Peak_EE = 0.15; Peak_I = 0.6;

% spatial indexes of E and I neurons
[NnS.X,NnS.Y] = V1Field_Generation(N_HC,1:N_S,'e',n_S_HC,n_I_HC);
[NnC.X,NnC.Y] = V1Field_Generation(N_HC,1:N_C,'e',n_C_HC,n_I_HC);
[NnI.X,NnI.Y] = V1Field_Generation(N_HC,1:N_I,'i');
%NnE.X = [NnS.X,NnC.X]; NnE.Y = [NnS.Y,NnC.Y];
%% 
% We are fixing boundary effect by connecting neurons to similar ODs

% Extended spatial indexes of E and I neurons (1HC each side)
Nex_S = n_S_HC^2 * (N_HC+2)^2; % neuron numbers In all
Nex_C = n_C_HC^2 * (N_HC+2)^2;
Nex_I = n_I_HC^2 * (N_HC+2)^2;
[NnSex.X,NnSex.Y] = V1Field_Generation(N_HC+2,1:Nex_S,'e',n_S_HC,n_I_HC);
[NnCex.X,NnCex.Y] = V1Field_Generation(N_HC+2,1:Nex_C,'e',n_C_HC,n_I_HC);
[NnIex.X,NnIex.Y] = V1Field_Generation(N_HC+2,1:Nex_I,'i');
%% 
% Now select complx cells

N_EE = 190; N_EI = 85; %*0.75
N_IE = 844; N_II = 84; %*0.75 Plan A:Original

CplxNEIncr = 255/N_EE; % Was 250 for 071021 simulation
EcplxInd = ismember(1:N_E, N_S+1:N_E);

N_SS = N_EE * (1-CplxR);            N_SC = N_EE * (CplxR);
N_CS = N_EE * (1-CplxR)*CplxNEIncr; N_CC = N_EE * (CplxR)*CplxNEIncr;
N_SI = N_EI;                        N_CI = N_EI;
N_IS = N_IE * (1-CplxR);            N_IC = N_IE * (CplxR);
%% Get a mapping between original and extended spatial indexs

if ~isfile([SaveFolder 'Paper2DriveNW_Conn.mat'])
    
    % From original to extended
    SMap_Ori2Ext = find(NnSex.X>n_S_HC & NnSex.X<=4*n_S_HC ...
        & NnSex.Y>n_S_HC & NnSex.Y<=4*n_S_HC);
    CMap_Ori2Ext = find(NnCex.X>n_C_HC & NnCex.X<=4*n_C_HC ...
        & NnCex.Y>n_C_HC & NnCex.Y<=4*n_C_HC);
    IMap_Ori2Ext = find(NnIex.X>n_I_HC & NnIex.X<=4*n_I_HC ...
        & NnIex.Y>n_I_HC & NnIex.Y<=4*n_I_HC);
    % from extended back to original: fold
    % E
    SMap_Ext2Ori = FoldExMap2Ori(NnSex,NnS,SMap_Ori2Ext,N_S,n_S_HC);
    CMap_Ext2Ori = FoldExMap2Ori(NnCex,NnC,CMap_Ori2Ext,N_C,n_C_HC);
    % I
    IMap_Ext2Ori = FoldExMap2Ori(NnIex,NnI,IMap_Ori2Ext,N_I,n_I_HC);
    % determine connections between E, I sparse metrices containing 0 or 1.
    % Row_i Column_j means neuron j projects to neuron i add periodic boundary
    % EE
    C_SS_Fix_Bd = ConnectionMat_FB_SimCplx(N_S,NnSex,Size_S, SMap_Ori2Ext, 1:N_S,...%
        N_S,NnSex,Size_S, SMap_Ext2Ori, 1:N_S,...
        Peak_EE,SD_E,Dist_LB,1, round(N_SS));
    C_SC_Fix_Bd = ConnectionMat_FB_SimCplx(N_S,NnSex,Size_S, SMap_Ori2Ext, 1:N_S,...%
        N_C,NnCex,Size_C, CMap_Ext2Ori, 1:N_C,...
        Peak_EE,SD_E,Dist_LB,0, round(N_SC));
    C_CS_Fix_Bd = ConnectionMat_FB_SimCplx(N_C,NnCex,Size_C, CMap_Ori2Ext, 1:N_C,...%
        N_S,NnSex,Size_S, SMap_Ext2Ori, 1:N_S,...
        Peak_EE,SD_E,Dist_LB,0, round(N_CS));
    C_CC_Fix_Bd = ConnectionMat_FB_SimCplx(N_C,NnCex,Size_C, CMap_Ori2Ext, 1:N_C,...%
        N_C,NnCex,Size_C, CMap_Ext2Ori, 1:N_C,...
        Peak_EE,SD_E,Dist_LB,1, round(N_CC));
    C_EE_Fix_Bd = [C_SS_Fix_Bd, C_SC_Fix_Bd; C_CS_Fix_Bd, C_CC_Fix_Bd];
    % EI
    C_SI_Fix_Bd = ConnectionMat_FB_SimCplx(N_S,NnSex,Size_S, SMap_Ori2Ext, 1:N_S,...
        N_I,NnIex,Size_I, IMap_Ext2Ori, 1:N_I,...
        Peak_I,SD_I,Dist_LB,0, round(N_SI));
    C_CI_Fix_Bd = ConnectionMat_FB_SimCplx(N_C,NnCex,Size_C, CMap_Ori2Ext, 1:N_C,...
        N_I,NnIex,Size_I, IMap_Ext2Ori, 1:N_I,...
        Peak_I,SD_I,Dist_LB,0, round(N_CI));
    C_EI_Fix_Bd = [C_SI_Fix_Bd; C_CI_Fix_Bd];
    % IE
    C_IS_Fix_Bd = ConnectionMat_FB_SimCplx(N_I,NnIex,Size_I, IMap_Ori2Ext, 1:N_I,...
        N_S,NnSex,Size_S, SMap_Ext2Ori, 1:N_S,...
        Peak_I,SD_E,Dist_LB,0, round(N_IS));
    C_IC_Fix_Bd = ConnectionMat_FB_SimCplx(N_I,NnIex,Size_I, IMap_Ori2Ext, 1:N_I,...
        N_C,NnCex,Size_C, CMap_Ext2Ori, 1:N_C,...
        Peak_I,SD_E,Dist_LB,0, round(N_IC));
    C_IE_Fix_Bd = [C_IS_Fix_Bd, C_IC_Fix_Bd];
    % II
    C_II_Fix_Bd = ConnectionMat_Fix_Boundary(N_I,NnIex,Size_I, IMap_Ori2Ext,...
        N_I,NnIex,Size_I, IMap_Ext2Ori,...
        Peak_I,SD_I,Dist_LB,1, round(N_II));
    
    CMatAll = struct('C_EE_Fix_Bd',C_EE_Fix_Bd,...
        'C_EI_Fix_Bd',C_EI_Fix_Bd,...
        'C_IE_Fix_Bd',C_IE_Fix_Bd,...
        'C_II_Fix_Bd',C_II_Fix_Bd);
else
    load([SaveFolder 'Paper2DriveNW_Conn.mat'])
    C_EE_Fix_Bd = CMatAll.C_EE_Fix_Bd ;
    C_EI_Fix_Bd = CMatAll.C_EI_Fix_Bd ;
    C_IE_Fix_Bd = CMatAll.C_IE_Fix_Bd ;
    C_II_Fix_Bd = CMatAll.C_II_Fix_Bd ;
end
%% 2. Variables and Parameters
% fixed parameters

DataPoint = 'V4D2'; %'V4D1'; 'V5D1'
load(sprintf("AllMFPixPara_Paper2TuneFig1%s.mat",DataPoint))

SaveFolder = [CurrentFolder '/Data/Paper2_NetworkTuning/Fig1V4/Paper2NWSimulationData/'];
close all
%% A: Now, with drive+LDE
% A1. Prepair parameters for drive simulation
% First figure out orientation domains

ODNum = 4; % 4 orientation domains
% A function from Neuron Ind and spatial scales to Orientation Domains
OD_S = zeros(size(NnS.X),'single');
OD_SMap = zeros(n_S_HC*N_HC,'single');
for NeuInd = 1:length(NnS.X)
    OD_S(NeuInd) = OrientDom(ODNum,NnS.X(NeuInd),NnS.Y(NeuInd),n_S_HC);
    OD_SMap(NnS.Y(NeuInd),NnS.X(NeuInd)) = OD_S(NeuInd);
end

OD_C = zeros(size(NnC.X),'single');
OD_CMap = zeros(n_C_HC*N_HC,'single');
for NeuInd = 1:length(NnC.X)
    OD_C(NeuInd) = OrientDom(ODNum,NnC.X(NeuInd),NnC.Y(NeuInd),n_C_HC);
    OD_CMap(NnC.Y(NeuInd),NnC.X(NeuInd)) = OD_C(NeuInd);
end
OD_E = [OD_S, OD_C];

OD_I = zeros(size(NnI.X),'single');
OD_IMap = zeros(n_I_HC*N_HC,'single');
for NeuInd = 1:length(NnI.X)
    OD_I(NeuInd) = OrientDom(ODNum,NnI.X(NeuInd),NnI.Y(NeuInd),n_I_HC);
    OD_IMap(NnI.Y(NeuInd),NnI.X(NeuInd)) = OD_I(NeuInd);
end
% Incorporate Phase:
% Use an updating vector for each E cell to reflect its on-off phase
% 
% Directly pull lgnEevents and do multiplication with the phase factor

StimulusFac = 1;
LGNL6Mapctgr = 3; % use cos maps
LGNFreq = 4; % 2 4 10Hz
tMod = 1e3/LGNFreq; MaxOnPhase = tMod/2;% Each Cycle lasts for 500 ms
PhaseE = single(tMod*rand(N_E,1));
%PhaseFRS = [90 67.5 45 67.5,    0  22.5 45 22.5] *N_Slgn/1e3 * StimulusFac; % 4 to 5
%        Both L6 and LGN: Should be blurred, and cos

Grating = 0; % degree
GratingHC = Grating + [0,45,90,135]; % grating for opt, obl1, ort, and obl2
L6up = 60; L6low = 6; %[10 54; 12 50]
% below: rates are already normalized to Hz/ms
if LGNL6Mapctgr == 1
    PhaseFRS_Grating = [45+abs(mod(GratingHC,180)-90)/2, 45-abs(mod(GratingHC,180)-90)/2]/1e3;
    L6Ord_F_Grating = ((abs(mod(GratingHC,180)-90)/90)      *(L6up-L6low)+L6low) /1e3; % get L6 frs
    ExpTex = 'Linear';
elseif LGNL6Mapctgr == 2
    PhaseFRS_Grating = [45+abs(mod(GratingHC,180)-90)/2, 45-abs(mod(GratingHC,180)-90)/2]/1e3;
    L6Ord_F_Grating = ((cosd(abs(mod(GratingHC,180))*2)+1)/2*(L6up-L6low)+L6low) /1e3;
    ExpTex = 'Cosine_L6';
elseif LGNL6Mapctgr == 3
    PhaseFRS_Grating = [45+(cosd(abs(mod(GratingHC,180))*2)+1)/2*45, 45-(cosd(abs(mod(GratingHC,180))*2)+1)/2*45]/1e3;
    L6Ord_F_Grating = ((cosd(abs(mod(GratingHC,180))*2)+1)/2*(L6up-L6low)+L6low) /1e3;
    
    ExpTex = 'Cosine_All';
else
    disp("illigal LGN->L6 mapping.")
end
%% 
% Filter L6 and LGN separately

% Map angles to LGN inpt
FigOn = false;
TruncLGN = 1.0; % truncation at X sigma
TruncL6 = 1.25;

sigLGN = 0.2; % was 0.1 for...
LGNFilt_Grating = SpatialGaussianFilt_my(OD_SMap,...
    N_HC,n_S_HC,n_S_HC*sigLGN,TruncLGN,FigOn);

LGNSon_Drive = LGNFilt_Grating * PhaseFRS_Grating(1:4)' *N_Slgn * StimulusFac;
LGNSoff_Drive = LGNFilt_Grating * PhaseFRS_Grating(5:8)' *N_Slgn * StimulusFac;

PhaseFRC = ones(1,8) *N_Clgn*(45/1e3) * StimulusFac;
lambda_EOn_drive  = [LGNSon_Drive; PhaseFRC(OD_C)']; % Should do another for Cplx cells
lambda_EOff_drive = [LGNSoff_Drive;PhaseFRC(OD_C+ODNum)'];
LGNE_Drive = [lambda_EOn_drive,lambda_EOff_drive];

% set a pre: E including both S and C?
lambda_E_drive_Pre = 45 *N_Slgn/1e3 * StimulusFac;
lambda_I_drive_Pre = 45 *N_Ilgn/1e3 * StimulusFac;

% Map angles to L6 inpt
sigL6 = 0.34;% of HC
L6Filt_GratingS = SpatialGaussianFilt_my(OD_SMap,...
    N_HC,n_S_HC,n_S_HC*sigL6,TruncL6,FigOn);
L6Filt_GratingC = SpatialGaussianFilt_my(OD_CMap,...
    N_HC,n_C_HC,n_C_HC*sigL6,TruncL6,FigOn);
L6Filt_GratingI = SpatialGaussianFilt_my(OD_IMap,...
    N_HC,n_I_HC,n_I_HC*sigL6,TruncL6,FigOn);

L6S_Drive = L6Filt_GratingS * L6Ord_F_Grating' *NS_L6 * StimulusFac;
L6C_Drive = L6Filt_GratingC * L6Ord_F_Grating' *NC_L6 * StimulusFac;
L6I_Drive = L6Filt_GratingI * L6Ord_F_Grating' *NI_L6 * StimulusFac;

rE_L6_Drive = [L6S_Drive;L6C_Drive];
rI_L6_Drive = L6I_Drive;

%% 
% LGN & L6 input classified by domains
% 
% Inputs are high now!!! Bernouli may underestimate, and we need to use:
%% 
% # first simulate Poisson processes
% # Feed them to the new function

% Making Poisson here/Reading
% T partition
T = 15000; dt = 0.1; TPar = 15;
TimeFrac = 0.05;
LGNCurInp = 0;

lgnE_Events = PoissonInputForNetwork(N_E,lambda_E_drive_Pre*(1-LGNCurInp),T*TimeFrac,dt);
lgnI_Events = PoissonInputForNetwork(N_I,lambda_I_drive_Pre*(1-LGNCurInp),T*TimeFrac,dt);
AmbE_Events = PoissonInputForNetwork(N_E,rE_amb,T*TimeFrac,dt);
AmbI_Events = PoissonInputForNetwork(N_I,rI_amb,T*TimeFrac,dt);

% L6 is much more specific!!! We will make a sum of Poisson and constant
L6CurInp = 0;
L6E_Events  = PoissonInputForNetwork(N_E,rE_L6_Drive*(1-L6CurInp),T*TimeFrac,dt);
L6I_Events  = PoissonInputForNetwork(N_I,rI_L6_Drive*(1-L6CurInp),T*TimeFrac,dt);
% A2. Simulation of Drive
% Then make a transfer between neurons in HCs and 10*10 pixels

%% Network Simulation
% SimulationT = 5000;
% E-to-E delay time
T_EEDly = .1; N_EEDly = floor(T_EEDly/dt);
T_IEDly = .1; N_IEDly = floor(T_IEDly/dt); %% NOTE: Was both 0.5
T_EIDly = .1; N_EIDly = floor(T_EIDly/dt);
% The initial states, for more simulation
load(['LargeNWFixIni.mat'],'EndState');
%load('IniTest.mat');
Fields = {'RefTimeE','VE','SpE','GE_ampa_R','GE_nmda_R','GE_gaba_R','GE_ampa_D','GE_nmda_D','GE_gaba_D',...
    'RefTimeI','VI','SpI','GI_ampa_R','GI_nmda_R','GI_gaba_R','GI_ampa_D','GI_nmda_D','GI_gaba_D'};
InSStr = cell2struct(EndState(1:length(Fields)), Fields, 2);
[InEs, InIs] = LargeNW_LoadIniState('InSStr', InSStr,N_E, N_I);
RefTimeE = InEs.RefTimeE; VE = InEs.VE; SpE = InEs.SpE; GE_ampa_R = InEs.GE_ampa_R; GE_nmda_R = InEs.GE_nmda_R; GE_gaba_R = InEs.GE_gaba_R;
GE_ampa_D = InEs.GE_ampa_D; GE_nmda_D = InEs.GE_nmda_D; GE_gaba_D = InEs.GE_gaba_D;
RefTimeI = InIs.RefTimeI; VI = InIs.VI; SpI = InIs.SpI; GI_ampa_R = InIs.GI_ampa_R; GI_nmda_R = InIs.GI_nmda_R; GI_gaba_R = InIs.GI_gaba_R;
GI_ampa_D = InIs.GI_ampa_D; GI_nmda_D = InIs.GI_nmda_D; GI_gaba_D = InIs.GI_gaba_D;
EEDlyRcd = sparse(N_E,N_EEDly);
IEDlyRcd = ones(N_I,N_IEDly);
EIDlyRcd = ones(N_E,N_EIDly);
% if N_EEDly>=size(EndState{19},2)
% EEDlyRcd = [EndState{19},zeros(N_E,N_EEDly-size(EndState{19},2))];
% else
%     EEDlyRcd = EndState{19}(:,1:N_EEDly);
% end
% if N_IEDly>=size(EndState{20},2)
% IEDlyRcd = [EndState{20},zeros(N_I,N_IEDly-size(EndState{20},2))];
% else
%     IEDlyRcd = EndState{20}(:,1:N_IEDly);
% end

% Create sliding windows
sampleT = 50;
Sliding = 50;
NSlide = floor(sampleT/Sliding);
TWinBounds = 0:Sliding:T;
Wins = [];
Wins(:,1) = TWinBounds(1:end-NSlide);
Wins(:,2) = TWinBounds(NSlide+1:end);
WinNum = length(Wins);
sampleN = floor(Sliding/dt); % sample each 2 ms
% SimulationN = floor(SimulationT/dt); % show and check every 200ms

BlowUp = false;
SampleInd = 1;

% Trace 1. Fr; 2. mV; 3. V distb
NWTrace = struct('SpEs',   cell(WinNum,1), 'SpIs',   cell(WinNum,1),...
    'mVEs',   cell(WinNum,1), 'mVIs',   cell(WinNum,1),...
    'GE_E',   cell(WinNum,1), 'GE_I',   cell(WinNum,1),...
    'GI_E',   cell(WinNum,1), 'GI_I',   cell(WinNum,1));
%NWTrace.Wins = {Wins};

MaxInpBin = floor(T*TimeFrac/dt);
%% 
% The biggest problem here is the access to the input matrices
% 
% So I do the following:
% 
% For every 100ms (1000 bins), first pre-compute input matrices for all (E ampa, 
% Enmda, Iampa, Inmda), then feed into iteration algorithm
% 
% Bug here: Though I defined AdjlgnE for every step, it is used only every 100ms...
% 
% So: 250on/250off are mixed here...

for TSec = 1:TPar
    tic
    for TimeN = 1:floor(T/dt/TPar)
        % Zero: Use Phase vec to determine different lgn phases for E neurons
        PhaseE = mod(PhaseE + dt,tMod);
        % First, Get input matrices from series
        InpWin = 100; FrameNum = floor(InpWin/dt);
        FrameInd = mod(TimeN, FrameNum);
        if FrameInd == 0
            FrameInd = FrameNum;
        end
        RandTimBin = randi([1 MaxInpBin+1],6,FrameNum);
        if FrameInd  == 1 % if the first frame, recompute input mats
            tt = 0:dt:InpWin-dt;
            % Judge the phase: if on, go first column, if off then second
            PhaseAll = floor(mod(PhaseE + tt, tMod)/MaxOnPhase)*N_E + repmat((1:N_E)',1,length(tt));
            AdjlgnE = LGNE_Drive(PhaseAll);
            
            EampaInp = single(full(S_Elgn * (lgnE_Events(:,RandTimBin(1,:)) .*AdjlgnE/lambda_E_drive_Pre + AdjlgnE*dt*LGNCurInp)...
                + S_amb  * AmbE_Events(:,RandTimBin(3,:)) ...
                + S_EL6  * ( L6E_Events(:,RandTimBin(5,:)) + rE_L6_Drive*dt*L6CurInp) * rhoE_ampa)); %  * rhoE_ampa
            IampaInp = single(full(S_Ilgn * (lgnI_Events(:,RandTimBin(2,:)) + lambda_I_drive_Pre*dt*LGNCurInp) ...
                + S_amb  * AmbI_Events(:,RandTimBin(4,:)) ...
                + S_IL6  * ( L6I_Events(:,RandTimBin(6,:)) + rI_L6_Drive*dt*L6CurInp) * rhoI_ampa)); %  * rhoI_ampa
            
            EnmdaInp = single(full(S_EL6  * ( L6E_Events(:,RandTimBin(5,:)) + rE_L6_Drive*dt*L6CurInp) * rhoE_nmda)); %
            InmdaInp = single(full(S_IL6  * ( L6I_Events(:,RandTimBin(6,:)) + rI_L6_Drive*dt*L6CurInp) * rhoI_nmda));
            
        end
        
        
        %      lgnEVec = zeros(N_E,1);lgnIVec = zeros(N_I,1);
        %      AmbEVec = zeros(N_E,1);AmbIVec = zeros(N_I,1);
        %      L6EVec = zeros(N_E,1); L6IVec = zeros(N_I,1);
        [oRefTimeE,oVE,oSpE,oGE_ampa_R,oGE_nmda_R,oGE_gaba_R,... % Output
            oGE_ampa_D,oGE_nmda_D,oGE_gaba_D,...
            oRefTimeI,oVI,oSpI,oGI_ampa_R,oGI_nmda_R,oGI_gaba_R,...
            oGI_ampa_D,oGI_nmda_D,oGI_gaba_D,...
            oEEDlyRcd,oIEDlyRcd, oEIDlyRcd] = ... % A updated N*T Mat recoding the time of kicks taking effect
            V1NetworkUpdate_Ver3_Drive_PoissonRead_EEnIEnEIDelay(RefTimeE,VE,SpE,GE_ampa_R,GE_nmda_R,GE_gaba_R,... % These are input kept updating
            GE_ampa_D,GE_nmda_D,GE_gaba_D,...
            RefTimeI,VI,SpI,GI_ampa_R,GI_nmda_R,GI_gaba_R,...
            GI_ampa_D,GI_nmda_D,GI_gaba_D,...
            EEDlyRcd,IEDlyRcd, EIDlyRcd,... % A N*T Mat recoding the time of kicks taking effect
            C_EE_Fix_Bd,C_EI_Fix_Bd,C_IE_Fix_Bd,C_II_Fix_Bd,...
            S_EE,S_EI,S_IE,S_II,...
            tau_ampa_R,tau_ampa_D,tau_nmda_R,tau_nmda_D,tau_gaba_R,tau_gaba_D,tau_ref,... % time unit is ms
            dt,p_EEFail,...
            gL_E,Ve,rhoE_ampa,rhoE_nmda,...
            gL_I,Vi,rhoI_ampa,rhoI_nmda,...
            EampaInp(:,FrameInd), IampaInp(:,FrameInd),...
            EnmdaInp(:,FrameInd), InmdaInp(:,FrameInd));
        % iteration
        RefTimeE = oRefTimeE; VE = oVE;SpE = oSpE;GE_ampa_R = oGE_ampa_R; GE_nmda_R = oGE_nmda_R; GE_gaba_R = oGE_gaba_R;
        GE_ampa_D = oGE_ampa_D; GE_nmda_D = oGE_nmda_D; GE_gaba_D = oGE_gaba_D;
        RefTimeI = oRefTimeI; VI = oVI;SpI = oSpI;GI_ampa_R = oGI_ampa_R; GI_nmda_R = oGI_nmda_R; GI_gaba_R = oGI_gaba_R;
        GI_ampa_D = oGI_ampa_D; GI_nmda_D = oGI_nmda_D; GI_gaba_D = oGI_gaba_D;
        EEDlyRcd = oEEDlyRcd; IEDlyRcd = oIEDlyRcd; EIDlyRcd = oEIDlyRcd;
        
        % Record Every Time time window
        RecordNum = floor(sampleT/dt);
        RecordInd = mod(TimeN, RecordNum);
        SampleRate = 0.1; % Can't be too small for multiple phases
        if RecordInd == 1
            clear E_Sp I_Sp mVETemp mVITemp
            E_Sp = [];
            I_Sp = [];
            mVETemp = zeros(N_E,floor(RecordNum*SampleRate),'single');
            mVITemp = zeros(N_I,floor(RecordNum*SampleRate),'single');
            GE_ETemp = zeros(N_E,floor(RecordNum*SampleRate),'single');
            GE_ITemp = zeros(N_E,floor(RecordNum*SampleRate),'single');
            GI_ETemp = zeros(N_I,floor(RecordNum*SampleRate),'single');
            GI_ITemp = zeros(N_I,floor(RecordNum*SampleRate),'single');
            mVRecordInd = 1;
        elseif RecordInd == 0
            NWTrace(SampleInd).mVEs = mVETemp;
            NWTrace(SampleInd).mVIs = mVITemp;
            
            NWTrace(SampleInd).SpEs = E_Sp;
            NWTrace(SampleInd).SpIs = I_Sp;
            
            NWTrace(SampleInd).GE_I = GE_ITemp;
            NWTrace(SampleInd).GE_E = GE_ETemp;
            NWTrace(SampleInd).GI_I = GI_ITemp;
            NWTrace(SampleInd).GI_E = GI_ETemp;
            SampleInd = SampleInd + 1
            %sum(isnan(oVE))/N_E
            FrESNow = sum(ismember(E_Sp(:,1), find(~EcplxInd)))/N_S/sampleT*1000
            FrECNow = sum(ismember(E_Sp(:,1), find(EcplxInd)))/N_C/sampleT*1000
            FrINow = size(I_Sp,1)/N_I/sampleT*1000
            toc
        end
        E_Sp = single([E_Sp;[find(oSpE),ones(size(find(oSpE)))*TimeN*dt]]);
        I_Sp = single([I_Sp;[find(oSpI),ones(size(find(oSpI)))*TimeN*dt]]);
        
        if mod(TimeN,floor(1/SampleRate)) == 5
            mVETemp(:,mVRecordInd) = oVE;
            mVITemp(:,mVRecordInd) = oVI;
            % Precompute GE GI
            GE_I = 1/(tau_gaba_D-tau_gaba_R) * (GE_gaba_D - GE_gaba_R); % S_EI is included in amplitude of GE_gaba
            GE_E = 1/(tau_ampa_D-tau_ampa_R) * (GE_ampa_D - GE_ampa_R) + ...
                1/(tau_nmda_D-tau_nmda_R) * (GE_nmda_D - GE_nmda_R); %
            
            GI_I = 1/(tau_gaba_D-tau_gaba_R) * (GI_gaba_D - GI_gaba_R); % S_EI is included in amplitude of GE_gaba
            GI_E = 1/(tau_ampa_D-tau_ampa_R) * (GI_ampa_D - GI_ampa_R) + ...
                1/(tau_nmda_D-tau_nmda_R) * (GI_nmda_D - GI_nmda_R); %
            GE_ITemp(:,mVRecordInd) = GE_I;
            GE_ETemp(:,mVRecordInd) = GE_E;
            GI_ITemp(:,mVRecordInd) = GI_I;
            GI_ETemp(:,mVRecordInd) = GI_E;
            %record loop goes forward
            mVRecordInd = mVRecordInd + 1;
        end
        
        if sum(isnan(oVE))>0.80*N_E
            BlowUp = true;
            disp('warning!: Network exploded')
            break
        end
        
        % the end of iteration
    end
    toc
    text = [sprintf('DriveWkSp_SCSepa_Cconst_%dHz_Deg%.1f_',LGNFreq,Grating) ...
        num2str(TSec) 's_NewSmear.mat'];
    EndState = {RefTimeE, VE, SpE, GE_ampa_R, GE_nmda_R, GE_gaba_R, GE_ampa_D, GE_nmda_D, GE_gaba_D,...
        RefTimeI, VI, SpI, GI_ampa_R, GI_nmda_R, GI_gaba_R, GI_ampa_D, GI_nmda_D, GI_gaba_D,...
        EEDlyRcd, IEDlyRcd};
    save([SaveFolder text],'NWTrace','EndState','PhaseE','-v7.3')
    
    clear NWTrace
    NWTrace = struct('SpEs',   cell(WinNum,1), 'SpIs',   cell(WinNum,1),...
        'mVEs',   cell(WinNum,1), 'mVIs',   cell(WinNum,1),...
        'GE_E',   cell(WinNum,1), 'GE_I',   cell(WinNum,1),...
        'GI_E',   cell(WinNum,1), 'GI_I',   cell(WinNum,1));
end
%save([SaveFolder,'Paper2DriveNW_Conn.mat'],'CMatAll','EcplxInd');
%%
DriveDataFolder = SaveFolder;% ...
% sprintf('SIEMt%.3f_SEIMt%.3f/',S_IE/S_II,S_EI/S_EE)];
%addpath(DriveDataFolder)
SpE = []; SpI = [];
mVE = []; mVI = []; % These are the total data
GE_E = []; GE_I = [];
GI_E = []; GI_I = [];
PhaseEAll = [];
for TSec = 1:TPar
    %text = [sprintf('DriveWkSp_SCSepa_Cconst_%dHz',LGNFreq) num2str(TSec) 's.mat'];
    text = [sprintf('DriveWkSp_SCSepa_Cconst_%dHz_Deg%.1f_',LGNFreq,Grating) ...
        num2str(TSec) 's_NewSmear.mat'];
    load([DriveDataFolder text],'NWTrace');
    
    SpECurrent = []; SpICurrent = [];
    mVECurrent = []; mVICurrent = [];
    GE_ECurrent = [];GE_ICurrent = [];
    GI_ECurrent = [];GI_ICurrent = [];
    PhaseECurrent = [];
    for WinInd = 1:length(NWTrace)
        SpECurrent = [SpECurrent; NWTrace(WinInd).SpEs];
        SpICurrent = [SpICurrent; NWTrace(WinInd).SpIs];
        mVECurrent = [mVECurrent, NWTrace(WinInd).mVEs];
        mVICurrent = [mVICurrent, NWTrace(WinInd).mVIs];
        
        %     GE_ECurrent = [GE_ECurrent, NWTrace(WinInd).GE_E];
        %     GE_ICurrent = [GE_ICurrent, NWTrace(WinInd).GE_I];
        %     GI_ECurrent = [GI_ECurrent, NWTrace(WinInd).GI_E];
        %     GI_ICurrent = [GI_ICurrent, NWTrace(WinInd).GI_I];
    end
    SpECurrent(:,2) = SpECurrent(:,2) + (TSec-1)*T/TPar;
    SpICurrent(:,2) = SpICurrent(:,2) + (TSec-1)*T/TPar;
    SpE = [SpE;SpECurrent];
    SpI = [SpI;SpICurrent];
    mVE = [mVE,mVECurrent];
    mVI = [mVI,mVICurrent];
    
    %     GE_E = [GE_E, GE_ECurrent];
    %     GE_I = [GE_I, GE_ICurrent];
    %     GI_E = [GI_E, GI_ECurrent];
    %     GI_I = [GI_I, GI_ICurrent];
    
end
load([DriveDataFolder text],'PhaseE');
load([DriveDataFolder 'Paper2DriveNW_Conn.mat'],'CMatAll','EcplxInd');
% Get conn mats
C_EE_Fix_Bd = CMatAll.C_EE_Fix_Bd ;
C_EI_Fix_Bd = CMatAll.C_EI_Fix_Bd ;
C_IE_Fix_Bd = CMatAll.C_IE_Fix_Bd ;
C_II_Fix_Bd = CMatAll.C_II_Fix_Bd ;

C_SS_Fix_Bd = C_EE_Fix_Bd(~EcplxInd,~EcplxInd);
C_SC_Fix_Bd = C_EE_Fix_Bd(~EcplxInd,EcplxInd);
C_CS_Fix_Bd = C_EE_Fix_Bd(EcplxInd,~EcplxInd);
C_CC_Fix_Bd = C_EE_Fix_Bd(EcplxInd,EcplxInd);
C_SI_Fix_Bd = C_EI_Fix_Bd(~EcplxInd,:);
C_CI_Fix_Bd = C_EI_Fix_Bd(EcplxInd,:);
C_IS_Fix_Bd = C_IE_Fix_Bd(:,~EcplxInd);
C_IC_Fix_Bd = C_IE_Fix_Bd(:,EcplxInd);

WinNum = size(mVE,2);
StatWin = linspace(0,T,WinNum+1);
StatWinSize = StatWin(2) - StatWin(1);
%% 
% Collect pixelwise Frs
%% 
% # Frs for every neuron

WinSize = 9000;
Window = [T-WinSize T];
EffSize = Window(2)-max(0,Window(1));

ESPinWin = SpE(:,2)>=Window(1) & SpE(:,2)<=Window(2);
SpEGood = SpE(ESPinWin,1);
[SpCN, SpCI] = groupcounts(SpEGood);
FrE_temp = zeros(N_E,1);
FrE_temp(SpCI) = SpCN*1e3/(EffSize);


ISPinWin = (SpI(:,2)>=Window(1) & SpI(:,2)<=Window(2));
SpIGood = SpI(ISPinWin,1);
[SpCN, SpCI] = groupcounts(SpIGood);
FrI_temp = zeros(N_I,1);
FrI_temp(SpCI) = SpCN*1e3/EffSize;
FrI_Map = reshape(FrI_temp,N_HC*n_I_HC,N_HC*n_I_HC);
%% 
% Pixelwise Frs: BEFORE inh suppresions!

[FrSPixMat,NnSPixel] = NeuVec2Pixel(FrE_temp(~EcplxInd),NnS,NPixX*N_HC,NPixY*N_HC);
[FrCPixMat,NnCPixel] = NeuVec2Pixel(FrE_temp(EcplxInd),NnC,NPixX*N_HC,NPixY*N_HC);
[FrIPixMat,NnIPixel] = NeuVec2Pixel(FrI_temp,NnI,NPixX*N_HC,NPixY*N_HC);
%% 
% Inhibition suppresion

%Thrsld = 70; Highist = [100,95];
%FrI_temp = InhKill(FrI_temp, Thrsld, Highist);
%FrI_temp = InhKill(FrI_temp, Thrsld, Highist);
%% 
% Pixelwise total L4E/I input
% 
% Total L4E: Sum of total L4E spikes received by 1S, 1C and 1I
% 
% Should put I supression here if...But not now.

L4EE = C_EE_Fix_Bd * FrE_temp;
L4EI = C_EI_Fix_Bd * FrI_temp;
L4IE = C_IE_Fix_Bd * FrE_temp;
L4II = C_II_Fix_Bd * FrI_temp;

[L4SEPix,~] = NeuVec2Pixel(L4EE(~EcplxInd),NnS,NPixX*N_HC,NPixY*N_HC);
[L4CEPix,~] = NeuVec2Pixel(L4EE(EcplxInd), NnC,NPixX*N_HC,NPixY*N_HC);
[L4IEPix,~] = NeuVec2Pixel(L4IE,           NnI,NPixX*N_HC,NPixY*N_HC);

[L4SIPix,~] = NeuVec2Pixel(L4EI(~EcplxInd),NnS,NPixX*N_HC,NPixY*N_HC);
[L4CIPix,~] = NeuVec2Pixel(L4EI(EcplxInd), NnC,NPixX*N_HC,NPixY*N_HC);
[L4IIPix,~] = NeuVec2Pixel(L4II,           NnI,NPixX*N_HC,NPixY*N_HC);
%% 
% Put Input/Coupling/Frs together

NWSmlt.FS = reshape(FrSPixMat,NPixX*N_HC*NPixY*N_HC,1);
NWSmlt.FC = reshape(FrCPixMat,NPixX*N_HC*NPixY*N_HC,1);
NWSmlt.FI = reshape(FrIPixMat,NPixX*N_HC*NPixY*N_HC,1);

NWSmlt.L4E = reshape((L4SEPix+L4CEPix)+L4IEPix, NPixX*N_HC*NPixY*N_HC,1);
NWSmlt.L4I = reshape(L4SIPix+L4CIPix+L4IIPix, NPixX*N_HC*NPixY*N_HC,1);

PixL6Ctgr = LGNIndSpat(L6Filt_GratingS,1:4,NnSPixel,N_HC,N_HC,NPixX,NPixY);
PixLGNCtgr = LGNIndSpat(LGNFilt_Grating,1:4,NnSPixel,N_HC,N_HC,NPixX,NPixY);

PixInptCtgrUse  = zeros(PixNum,4,4);
PixNum = N_HC^2*NPixX*NPixY;
for PixInd = 1:PixNum
    PixInptCtgrUse(PixInd,:,:) = PixLGNCtgr(PixInd,:)' *PixL6Ctgr(PixInd,:);% by multiplying both indexes
end

NWSmlt.PixInptCtgrUse = PixInptCtgrUse;

save([SaveFolder sprintf('NWSimulationPix_%.1fdeg_NewSmear.mat',Grating)],'NWSmlt')
