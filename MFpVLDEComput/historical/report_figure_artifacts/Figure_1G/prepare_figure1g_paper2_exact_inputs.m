function prepare_figure1g_paper2_exact_inputs(bundleFile,paramFile,outputRoot)
% Recreate the external input files expected by the unmodified Paper 2 driver.
if ~exist(outputRoot,'dir'), mkdir(outputRoot); end
B=load(bundleFile,'C_EE_Fix_Bd','C_EI_Fix_Bd','C_IE_Fix_Bd','C_II_Fix_Bd');
P=load(paramFile,'EndState','EcplxInd');
CMatAll=struct('C_EE_Fix_Bd',B.C_EE_Fix_Bd, ...
    'C_EI_Fix_Bd',B.C_EI_Fix_Bd, ...
    'C_IE_Fix_Bd',B.C_IE_Fix_Bd, ...
    'C_II_Fix_Bd',B.C_II_Fix_Bd);
EcplxInd=P.EcplxInd;
EndState=P.EndState;
save(fullfile(outputRoot,'Paper2DriveNW_Conn.mat'),'CMatAll','EcplxInd','-v7.3');
save(fullfile(outputRoot,'LargeNWFixIni.mat'),'EndState','-v7.3');
end
