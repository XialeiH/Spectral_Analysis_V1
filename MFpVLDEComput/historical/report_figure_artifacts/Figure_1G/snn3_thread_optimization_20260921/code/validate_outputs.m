function validate_outputs(candidate,control,output)
% Outside official timers: exact final rates and all saved endpoint states.
suffix='Data/Paper2_NetworkTuning/Fig1V4/Paper2NWSimulationData';
a=fullfile(candidate,suffix); b=fullfile(control,suffix);
A=load(fullfile(a,'NWSimulationPix_0.0deg_NewSmear.mat'),'NWSmlt');
B=load(fullfile(b,'NWSimulationPix_0.0deg_NewSmear.mat'),'NWSmlt');
finalEqual=isequaln(A.NWSmlt,B.NWSmlt);
rateMax=max(abs([A.NWSmlt.FS-B.NWSmlt.FS;A.NWSmlt.FC-B.NWSmlt.FC;A.NWSmlt.FI-B.NWSmlt.FI]));
segmentEqual=false(20,1);
for k=1:20
    name=sprintf('DriveWkSp_SCSepa_Cconst_4Hz_Deg0.0_%ds_NewSmear.mat',k);
    A=load(fullfile(a,name),'EndState','PhaseE');
    B=load(fullfile(b,name),'EndState','PhaseE');
    segmentEqual(k)=isequaln(A,B);
end
result=struct('final_equal',finalEqual,'rate_max_absolute',rateMax, ...
    'segment_endpoints_equal',segmentEqual,'passed',finalEqual && all(segmentEqual));
f=fopen(output,'w'); fprintf(f,'%s\n',jsonencode(result)); fclose(f);
assert(result.passed,'Thread change alters SNN output; do not publish as equivalent.');
end
