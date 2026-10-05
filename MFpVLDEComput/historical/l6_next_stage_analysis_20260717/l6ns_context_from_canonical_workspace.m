% Assemble ModelContext after the canonical h96 one-condition setup script.
% This is a script intentionally: it reads the canonical variables from the
% caller workspace and leaves ModelContext there.

ModelContext = struct();
ModelContext.PixLGNCtgr = PixLGNCtgr;
ModelContext.L6Kernel = L6Kernel;
ModelContext.L6Parameters = L6pars{L6parId};
ModelContext.C_SS = C_SS_meanU;
ModelContext.C_CS = C_CS_meanU;
ModelContext.C_IS = C_IS_mean;
ModelContext.C_SC = C_SC_meanU;
ModelContext.C_CC = C_CC_meanU;
ModelContext.C_IC = C_IC_mean;
ModelContext.C_SI = C_SI_mean;
ModelContext.C_CI = C_CI_mean;
ModelContext.C_II = C_II_mean;
ModelContext.L4SEp = L4SEp;
ModelContext.L4SIp = L4SIp;
ModelContext.L4CEp = L4CEp;
ModelContext.L4CIp = L4CIp;
ModelContext.L4IEp = L4IEp;
ModelContext.L4IIp = L4IIp;
ModelContext.ContrastUse = Contruse;
ModelContext.OrientationUse = Orientationuse;
ModelContext.EKpUse = EKpUse;
ModelContext.IKpUse = IKpUse;
ModelContext.Isaturation = Isaturation;
ModelContext.CWeight = 0.3077;
ModelContext.RelaxationP = p;
ModelContext.InitialStateStruct = struct('S',IniTest.S(:), ...
    'C',IniTest.C(:),'I',IniTest.I(:));
ModelContext.InitialState = [IniTest.S(:);IniTest.C(:);IniTest.I(:)];
ModelContext.FixedPointStruct = struct('S',LDEfixedpoint.S(:), ...
    'C',LDEfixedpoint.C(:),'I',LDEfixedpoint.I(:));
ModelContext.FixedPoint = [LDEfixedpoint.S(:);LDEfixedpoint.C(:);LDEfixedpoint.I(:)];
ModelContext.MapSize = [N_HCOutY*NPixY, numel(LDEfixedpoint.S)/(N_HCOutY*NPixY)];

ModelContext.FixedL6LibInd = l6ns_l6_dynamic(ModelContext.FixedPoint,ModelContext,[1 1]);
