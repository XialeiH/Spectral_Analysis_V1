function run_followup_threshold_task()
% Test whether the branch jump follows changes to the extension threshold.

paths=followup_initialize();
loaded=load(paths.SetupFile,'setup');
context=loaded.setup.Context;
baseline=context.FixedPoint(:);
specs=followup_specs('threshold');
taskId=local_task_id();
spec=specs(taskId,:);
if spec.pathway=="L6"
    betaGrid=(0.22:0.005:0.35)';
else
    betaGrid=(-0.05:-0.0025:-0.14)';
end
continuationState=baseline;
rows=cell(2*numel(betaGrid),1);
for index=1:numel(betaGrid)
    beta=betaGrid(index);
    [gain6,gainI]=local_gains(spec.pathway,beta);
    if spec.mapMode=="raw"
        phi=@(x)followup_phi_threshold(x,context,gain6,gainI,"raw",1,1);
    else
        phi=@(x)followup_phi_threshold(x,context,gain6,gainI,"extended", ...
            spec.edgeScale,1);
    end
    resetResult=real_tuning_fixed_point(phi,baseline,context.RelaxationP);
    continuationResult=real_tuning_fixed_point(phi,continuationState, ...
        context.RelaxationP);
    if continuationResult.Converged
        continuationState=continuationResult.State(:);
    end
    rows{2*index-1}=local_row(taskId,spec,beta,"baseline reset", ...
        resetResult,phi,context,gain6,gainI);
    rows{2*index}=local_row(taskId,spec,beta,"low continuation", ...
        continuationResult,phi,context,gain6,gainI);
end
resultTable=vertcat(rows{:});
if spec.mapMode=="extended" && abs(spec.edgeScale-1)<1e-12
    testState=resultTable.maximumRateHz(end)>180;
    stateToTest=continuationState;
    gain6=1; gainI=1;
    if spec.pathway=="L6"; gain6=1+betaGrid(end); else; gainI=1+betaGrid(end); end
    exact=mechanism_phi_variant(stateToTest,context,gain6,gainI,'extended');
    custom=followup_phi_threshold(stateToTest,context,gain6,gainI,"extended",1,1);
    matchError=norm(exact-custom)/max(norm(exact),eps);
    fprintf('scale-one map match %.3e, high=%d\n',matchError,testState);
    if matchError>2e-10
        error('Followup:ThresholdMapMismatch','Scale-one map mismatch %.3e.',matchError);
    end
end
outDir=fullfile(paths.FollowupOutput,'threshold_tasks');
if ~exist(outDir,'dir'); mkdir(outDir); end
writetable(resultTable,fullfile(outDir,sprintf('threshold_%03d.tsv',taskId)), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outDir,sprintf('threshold_%03d.mat',taskId)), ...
    'resultTable','-v7.3');
fprintf('%s mode=%s edgeScale=%g: reset high at beta %s; continuation high at beta %s\n', ...
    spec.pathway,spec.mapMode,spec.edgeScale, ...
    local_first_high(resultTable,"baseline reset"), ...
    local_first_high(resultTable,"low continuation"));
end

function row=local_row(taskId,spec,beta,protocol,result,phi,context,gain6,gainI)
state=result.State(:);
residual=norm(phi(state)-state)/max(norm(state),eps);
[outside,meanWeight,fullFraction]=local_extension_stats( ...
    state,context,gain6,gainI,spec);
raw=mechanism_phi_variant(state,context,gain6,gainI,'raw');
rawResidual=norm(raw-state)/max(norm(state),eps);
row=table(taskId,spec.pathway,spec.mapMode,spec.edgeScale,beta,string(protocol), ...
    result.Converged,string(result.Termination),result.Iterations,residual, ...
    rawResidual,min(state),max(state),outside,meanWeight,fullFraction, ...
    max(state)>180, ...
    'VariableNames',{'taskId','pathway','mapMode','edgeScale','beta','protocol', ...
    'converged','termination','iterations','mapResidual','rawMapResidual', ...
    'minimumRateHz','maximumRateHz','fractionOutside','meanExtensionWeight', ...
    'fractionFullExtension','highBranch'});
end

function [outside,meanWeight,fullFraction]=local_extension_stats( ...
        state,context,gain6,gainI,spec)
if spec.mapMode=="raw"
    outside=0; meanWeight=0; fullFraction=0; return
end
inputs=mechanism_state_inputs(state,context,gain6,gainI);
edgeE=47500*spec.edgeScale; edgeI=35315.63671875*spec.edgeScale;
weights=[];
for name={'S','C','I'}
    E=inputs.(name{1}).L4E; I=inputs.(name{1}).L4I;
    distance=sqrt((max(0,E-edgeE)/5000).^2+(max(0,I-edgeI)/5000).^2);
    w=ones(size(distance)); mask=distance>0 & distance<1; t=distance(mask);
    w(mask)=6*t.^5-15*t.^4+10*t.^3; w(distance<=0)=0;
    weights=[weights;w]; %#ok<AGROW>
end
outside=mean(weights>0); meanWeight=mean(weights); fullFraction=mean(weights>=1);
end

function textValue=local_first_high(tableData,protocol)
rows=tableData(tableData.protocol==protocol & tableData.highBranch,:);
if isempty(rows); textValue='none'; else; textValue=sprintf('%+.4f',rows.beta(1)); end
end

function taskId=local_task_id()
taskId=str2double(getenv('FOLLOWUP_TASK_ID'));
if ~isfinite(taskId); taskId=str2double(getenv('SLURM_ARRAY_TASK_ID')); end
if ~isfinite(taskId); taskId=1; end
taskId=round(taskId);
end

function [gain6,gainI]=local_gains(pathway,beta)
gain6=1; gainI=1;
if pathway=="L6"; gain6=1+beta; else; gainI=1+beta; end
end
