function build_figure1g_snn3_bundle(repoRoot, outputFile)
% Prepare the fixed 3x3 Paper 2 SNN state outside benchmark timers.
addpath(fullfile(repoRoot, 'Utils'));
paramFile = fullfile(repoRoot, 'Data', 'Paper2_NetworkTuning', ...
    'Fig1V4', 'AllMFPixPara_Paper2TuneFig1V4D2.mat');
P = load(paramFile);
assert(P.N_HC == 3, 'The authoritative Paper 2 workspace is not 3x3.');
rng(20260901, 'twister');

fprintf('Building fixed 3x3 Paper 2 connectivity.\n');
C_SS = ConnectionMat_FB_SimCplx(P.N_S,P.NnSex,P.Size_S,P.SMap_Ori2Ext,1:P.N_S, ...
    P.N_S,P.NnSex,P.Size_S,P.SMap_Ext2Ori,1:P.N_S, ...
    P.Peak_EE,P.SD_E,P.Dist_LB,1,round(P.N_SS));
C_SC = ConnectionMat_FB_SimCplx(P.N_S,P.NnSex,P.Size_S,P.SMap_Ori2Ext,1:P.N_S, ...
    P.N_C,P.NnCex,P.Size_C,P.CMap_Ext2Ori,1:P.N_C, ...
    P.Peak_EE,P.SD_E,P.Dist_LB,0,round(P.N_SC));
C_CS = ConnectionMat_FB_SimCplx(P.N_C,P.NnCex,P.Size_C,P.CMap_Ori2Ext,1:P.N_C, ...
    P.N_S,P.NnSex,P.Size_S,P.SMap_Ext2Ori,1:P.N_S, ...
    P.Peak_EE,P.SD_E,P.Dist_LB,0,round(P.N_CS));
C_CC = ConnectionMat_FB_SimCplx(P.N_C,P.NnCex,P.Size_C,P.CMap_Ori2Ext,1:P.N_C, ...
    P.N_C,P.NnCex,P.Size_C,P.CMap_Ext2Ori,1:P.N_C, ...
    P.Peak_EE,P.SD_E,P.Dist_LB,1,round(P.N_CC));
C_EE_Fix_Bd = [C_SS,C_SC; C_CS,C_CC];
clear C_SS C_SC C_CS C_CC

C_SI = ConnectionMat_FB_SimCplx(P.N_S,P.NnSex,P.Size_S,P.SMap_Ori2Ext,1:P.N_S, ...
    P.N_I,P.NnIex,P.Size_I,P.IMap_Ext2Ori,1:P.N_I, ...
    P.Peak_I,P.SD_I,P.Dist_LB,0,round(P.N_SI));
C_CI = ConnectionMat_FB_SimCplx(P.N_C,P.NnCex,P.Size_C,P.CMap_Ori2Ext,1:P.N_C, ...
    P.N_I,P.NnIex,P.Size_I,P.IMap_Ext2Ori,1:P.N_I, ...
    P.Peak_I,P.SD_I,P.Dist_LB,0,round(P.N_CI));
C_EI_Fix_Bd = [C_SI; C_CI];
clear C_SI C_CI

C_IS = ConnectionMat_FB_SimCplx(P.N_I,P.NnIex,P.Size_I,P.IMap_Ori2Ext,1:P.N_I, ...
    P.N_S,P.NnSex,P.Size_S,P.SMap_Ext2Ori,1:P.N_S, ...
    P.Peak_I,P.SD_E,P.Dist_LB,0,round(P.N_IS));
C_IC = ConnectionMat_FB_SimCplx(P.N_I,P.NnIex,P.Size_I,P.IMap_Ori2Ext,1:P.N_I, ...
    P.N_C,P.NnCex,P.Size_C,P.CMap_Ext2Ori,1:P.N_C, ...
    P.Peak_I,P.SD_E,P.Dist_LB,0,round(P.N_IC));
C_IE_Fix_Bd = [C_IS,C_IC];
clear C_IS C_IC

C_II_Fix_Bd = ConnectionMat_Fix_Boundary(P.N_I,P.NnIex,P.Size_I,P.IMap_Ori2Ext, ...
    P.N_I,P.NnIex,P.Size_I,P.IMap_Ext2Ori, ...
    P.Peak_I,P.SD_I,P.Dist_LB,1,round(P.N_II));

StimulusFac = 1;
ODNum = 4;
LGNFreq = 4;
tMod = 1e3/LGNFreq;
MaxOnPhase = tMod/2;
GratingHC = [0,45,90,135];
L6up = 60;
L6low = 6;
PhaseFRS = [45+(cosd(abs(mod(GratingHC,180))*2)+1)/2*45, ...
    45-(cosd(abs(mod(GratingHC,180))*2)+1)/2*45]/1e3;
L6OrdF = ((cosd(abs(mod(GratingHC,180))*2)+1)/2*(L6up-L6low)+L6low)/1e3;
LGNFilt = SpatialGaussianFilt_my(P.OD_SMap,P.N_HC,P.n_S_HC, ...
    P.n_S_HC*0.2,1.0,false);
LGNSon = LGNFilt*PhaseFRS(1:4)'*P.N_Slgn*StimulusFac;
LGNSoff = LGNFilt*PhaseFRS(5:8)'*P.N_Slgn*StimulusFac;
PhaseFRC = ones(1,8)*P.N_Clgn*(45/1e3)*StimulusFac;
lambda_EOn = [LGNSon; PhaseFRC(P.OD_C)'];
lambda_EOff = [LGNSoff; PhaseFRC(P.OD_C+ODNum)'];
LGNE_Drive = [lambda_EOn,lambda_EOff];

L6FiltS = SpatialGaussianFilt_my(P.OD_SMap,P.N_HC,P.n_S_HC, ...
    P.n_S_HC*0.34,1.25,false);
L6FiltC = SpatialGaussianFilt_my(P.OD_CMap,P.N_HC,P.n_C_HC, ...
    P.n_C_HC*0.34,1.25,false);
L6FiltI = SpatialGaussianFilt_my(P.OD_IMap,P.N_HC,P.n_I_HC, ...
    P.n_I_HC*0.34,1.25,false);
rE_L6_Drive = [L6FiltS*L6OrdF'*P.NS_L6; L6FiltC*L6OrdF'*P.NC_L6];
rI_L6_Drive = L6FiltI*L6OrdF'*P.NI_L6;

Fields = {'RefTimeE','VE','SpE','GE_ampa_R','GE_nmda_R','GE_gaba_R', ...
    'GE_ampa_D','GE_nmda_D','GE_gaba_D','RefTimeI','VI','SpI', ...
    'GI_ampa_R','GI_nmda_R','GI_gaba_R','GI_ampa_D','GI_nmda_D','GI_gaba_D'};
InSStr = cell2struct(P.EndState(1:numel(Fields)),Fields,2);
[InEs,InIs] = LargeNW_LoadIniState('InSStr',InSStr,P.N_E,P.N_I);

B = struct;
B.N_HC = P.N_HC;
B.N_E = P.N_E;
B.N_I = P.N_I;
B.C_EE_Fix_Bd = C_EE_Fix_Bd;
B.C_EI_Fix_Bd = C_EI_Fix_Bd;
B.C_IE_Fix_Bd = C_IE_Fix_Bd;
B.C_II_Fix_Bd = C_II_Fix_Bd;
B.InEs = InEs;
B.InIs = InIs;
B.LGNE_Drive = LGNE_Drive;
B.rE_L6_Drive = rE_L6_Drive;
B.rI_L6_Drive = rI_L6_Drive;
B.tMod = tMod;
B.MaxOnPhase = MaxOnPhase;
names = {'S_EE','S_EI','S_IE','S_II','tau_ampa_R','tau_ampa_D', ...
    'tau_nmda_R','tau_nmda_D','tau_gaba_R','tau_gaba_D','tau_ref', ...
    'p_EEFail','gL_E','Ve','rhoE_ampa','rhoE_nmda','gL_I','Vi', ...
    'rhoI_ampa','rhoI_nmda','S_Elgn','S_Ilgn','S_amb','S_EL6','S_IL6', ...
    'N_Slgn','N_Ilgn','rE_amb','rI_amb'};
for k = 1:numel(names)
    B.(names{k}) = P.(names{k});
end
B.lambda_E_drive_Pre = 45*P.N_Slgn/1e3*StimulusFac;
B.lambda_I_drive_Pre = 45*P.N_Ilgn/1e3*StimulusFac;
B.dt = 0.1;
B.T = 10000;
B.TPar = 10;
B.TimeFrac = 0.05;
B.LGNCurInp = 0;
B.L6CurInp = 0;
B.sourceCommit = '9338749de32a3ab10340430ff1604c7135064b65';
B.sourceDriver = 'Paper2_Fig7Comp_NW_LDE.mlx';
B.bundleSeed = 20260901;
save(outputFile,'-struct','B','-v7.3');
fprintf('Saved %s\n',outputFile);
end
