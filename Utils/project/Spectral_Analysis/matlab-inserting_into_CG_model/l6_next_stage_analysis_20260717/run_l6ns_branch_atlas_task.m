function run_l6ns_branch_atlas_task
% Continue one catalog seed and direction selected by SLURM_ARRAY_TASK_ID.

root=getenv('L6NS_ATLAS_ROOT');
if isempty(root); root=pwd; end
taskId=str2double(getenv('SLURM_ARRAY_TASK_ID'));
if ~isfinite(taskId); taskId=1; end

loaded=load(fullfile(root,'results_global_bifurcation_20260722', ...
    'global_bifurcation_setup.mat'));
names=fieldnames(loaded);
setup=loaded.(names{1});
loaded=load(fullfile(root,'branch_atlas_20260803','atlas_seed_catalog.mat'), ...
    'catalog');
catalog=loaded.catalog;
eligibleClasses=["census_root","branch_intersection", ...
    "nearsilent_family","nearsilent_fold_seed"];
eligible=find(ismember(catalog.Table.seedClass,eligibleClasses));
numberOfSeeds=numel(eligible);
if taskId<1 || taskId>2*numberOfSeeds
    error('Task %d is outside the valid range 1:%d.',taskId,2*numberOfSeeds);
end
seedPosition=ceil(taskId/2);
direction=2*mod(taskId,2)-1;
seedIndex=eligible(seedPosition);
row=catalog.Table(seedIndex,:);
outputDir=fullfile(root,'branch_atlas_20260803','continuation_tasks');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
outputFile=fullfile(outputDir,sprintf('atlas_task_%03d_seed_%03d_dir_%+d.mat', ...
    taskId,row.seedId,direction));
maximumPoints=local_environment_number('L6NS_MAX_POINTS',140);
modeCount=local_environment_number('L6NS_MODE_COUNT',32);
options=struct('MinimumStep',1e-4,'MaximumStep',0.06, ...
    'InitialStep',0.008,'MaximumPoints',maximumPoints, ...
    'WeightBounds',[-1 2],'EtaW',10,'ModeCount',modeCount);
result=l6ns_atlas_weighted_continue(setup,char(row.pathway), ...
    catalog.States{seedIndex},row.freezeWeight,direction,outputFile,options);
metadata=table(taskId,row.seedId,row.pathway,row.freezeWeight,row.label, ...
    row.seedClass,direction,height(result.Summary),string(result.TerminalReason), ...
    'VariableNames',{'taskId','seedId','pathway','seedWeight','seedLabel', ...
    'seedClass','direction','pointCount','terminalReason'});
writetable(metadata,fullfile(outputDir,sprintf('atlas_task_%03d_metadata.tsv',taskId)), ...
    'FileType','text','Delimiter','\t');
end

function value=local_environment_number(name,defaultValue)
textValue=getenv(name);
value=str2double(textValue);
if isempty(textValue) || ~isfinite(value) || value<=0
    value=defaultValue;
end
end
