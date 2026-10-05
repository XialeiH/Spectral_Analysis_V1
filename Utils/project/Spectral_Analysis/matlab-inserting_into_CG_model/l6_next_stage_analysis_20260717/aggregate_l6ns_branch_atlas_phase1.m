function atlas = aggregate_l6ns_branch_atlas_phase1(root)
% Reconcile continuation tasks into connected components and event candidates.

if nargin<1 || isempty(root); root=pwd; end
taskDir=fullfile(root,'branch_atlas_20260803','continuation_tasks');
files=dir(fullfile(taskDir,'atlas_task_*_seed_*_dir_*.mat'));
if isempty(files); error('No continuation task files found in %s.',taskDir); end
tasks=cell(numel(files),1);
taskRows={};
pointTables=cell(numel(files),1);
eventRows={};
for fileIndex=1:numel(files)
    loaded=load(fullfile(files(fileIndex).folder,files(fileIndex).name),'result');
    task=loaded.result;
    % Component reconciliation uses summaries and fixed-point states only.
    % Do not retain per-point mode tables/vectors from all tasks in memory.
    if isfield(task,'ModeTables'); task=rmfield(task,'ModeTables'); end
    if isfield(task,'ModeVectors'); task=rmfield(task,'ModeVectors'); end
    tasks{fileIndex}=task;
    token=regexp(files(fileIndex).name, ...
        'atlas_task_(\d+)_seed_(\d+)_dir_([+-]\d+)','tokens','once');
    taskId=str2double(token{1}); seedId=str2double(token{2});
    direction=str2double(token{3});
    taskRows(end+1,:)={fileIndex,taskId,seedId,string(task.Pathway), ... %#ok<AGROW>
        direction,height(task.Summary),string(task.TerminalReason),files(fileIndex).name};
    points=task.Summary;
    points.taskIndex=repmat(fileIndex,height(points),1);
    points.taskId=repmat(taskId,height(points),1);
    points.seedId=repmat(seedId,height(points),1);
    points.direction=repmat(direction,height(points),1);
    points.branchPoint=(1:height(points))';
    pointTables{fileIndex}=points;
    for pointIndex=1:height(points)-1
        first=points(pointIndex,:); second=points(pointIndex+1,:);
        if first.dwds*second.dwds<0
            eventRows(end+1,:)=local_event_row(fileIndex,taskId,seedId, ... %#ok<AGROW>
                task.Pathway,'tangent_reversal',pointIndex,first,second);
        end
        if first.unstableDimension~=second.unstableDimension
            eventRows(end+1,:)=local_event_row(fileIndex,taskId,seedId, ... %#ok<AGROW>
                task.Pathway,'morse_index_change',pointIndex,first,second);
        elseif (first.maxRealLambda-1)*(second.maxRealLambda-1)<0
            eventRows(end+1,:)=local_event_row(fileIndex,taskId,seedId, ... %#ok<AGROW>
                task.Pathway,'stability_boundary_crossing',pointIndex,first,second);
        end
    end
end
taskTable=cell2table(taskRows,'VariableNames',{'taskIndex','taskId','seedId', ...
    'pathway','direction','pointCount','terminalReason','fileName'});
points=vertcat(pointTables{:});
if isempty(eventRows)
    events=table();
else
    events=cell2table(eventRows,'VariableNames',{'taskIndex','taskId','seedId', ...
        'pathway','eventType','leftPoint','rightPoint','leftWeight','rightWeight', ...
        'leftMaxRealLambda','rightMaxRealLambda','leftUnstableDimension', ...
        'rightUnstableDimension','leftDwds','rightDwds','minimumSigmaMinA'});
end

adjacency=speye(numel(tasks));
for firstTask=1:numel(tasks)-1
    for secondTask=firstTask+1:numel(tasks)
        if taskTable.pathway(firstTask)~=taskTable.pathway(secondTask); continue; end
        if taskTable.seedId(firstTask)==taskTable.seedId(secondTask)
            adjacency(firstTask,secondTask)=1; adjacency(secondTask,firstTask)=1;
            continue
        end
        if local_tasks_intersect(tasks{firstTask},tasks{secondTask})
            adjacency(firstTask,secondTask)=1; adjacency(secondTask,firstTask)=1;
        end
    end
end
componentId=conncomp(graph(adjacency),'Type','weak')';
taskTable.componentId=componentId;
points.componentId=componentId(points.taskIndex);
componentRows={};
for component=unique(componentId)'
    indices=find(componentId==component);
    componentPoints=points(ismember(points.taskIndex,indices),:);
    stablePoints=componentPoints(componentPoints.maxRealLambda<1,:);
    stableHighRatePoints=stablePoints(stablePoints.extremeHighRate,:);
    if isempty(stablePoints)
        maximumStableRateRatio=NaN;
        maximumStableRate=NaN;
    else
        maximumStableRateRatio=max(stablePoints.maximumRateRatio);
        maximumStableRate=max(stablePoints.maximumState);
    end
    componentRows(end+1,:)={component,taskTable.pathway(indices(1)), ... %#ok<AGROW>
        numel(indices),numel(unique(taskTable.seedId(indices))), ...
        min(componentPoints.freezeWeight),max(componentPoints.freezeWeight), ...
        min(componentPoints.maxRealLambda),max(componentPoints.maxRealLambda), ...
        min(componentPoints.unstableDimension),max(componentPoints.unstableDimension), ...
        min(componentPoints.sigmaMinA),height(stablePoints), ...
        height(stableHighRatePoints),maximumStableRateRatio,maximumStableRate, ...
        ~isempty(stableHighRatePoints)};
end
components=cell2table(componentRows,'VariableNames',{'componentId','pathway', ...
    'taskCount','seedCount','minimumWeight','maximumWeight','minimumMaxRealLambda', ...
    'maximumMaxRealLambda','minimumUnstableDimension','maximumUnstableDimension', ...
    'minimumSigmaMinA','stablePointCount','stableHighRatePointCount', ...
    'maximumStableRateRatio','maximumStableRate','containsStableHighRateState'});

atlas=struct('Tasks',taskTable,'Points',points,'Events',events, ...
    'Components',components,'Adjacency',adjacency,'TaskResults',{tasks});
outputDir=fullfile(root,'branch_atlas_20260803','aggregate');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
save(fullfile(outputDir,'phase1_branch_atlas.mat'),'atlas','-v7.3');
writetable(taskTable,fullfile(outputDir,'continuation_tasks.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(points,fullfile(outputDir,'continuation_points.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(events,fullfile(outputDir,'event_candidates.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(components,fullfile(outputDir,'connected_components.tsv'), ...
    'FileType','text','Delimiter','\t');
end

function row=local_event_row(taskIndex,taskId,seedId,pathway,eventType,index,first,second)
row={taskIndex,taskId,seedId,string(pathway),string(eventType),index,index+1, ...
    first.freezeWeight,second.freezeWeight,first.maxRealLambda, ...
    second.maxRealLambda,first.unstableDimension,second.unstableDimension, ...
    first.dwds,second.dwds,min(first.sigmaMinA,second.sigmaMinA)};
end

function intersects=local_tasks_intersect(first,second)
intersects=false;
a=first.Summary; b=second.Summary;
featureA=[10*a.freezeWeight,a.meanS/10,a.meanC/50,a.meanI/80, ...
    a.rmsDistanceS/30,a.rmsDistanceC/60,a.rmsDistanceI/80];
featureB=[10*b.freezeWeight,b.meanS/10,b.meanC/50,b.meanI/80, ...
    b.rmsDistanceS/30,b.rmsDistanceC/60,b.rmsDistanceI/80];
distanceSquared=max(0,sum(featureA.^2,2)+sum(featureB.^2,2)'-2*(featureA*featureB'));
[values,linearIndices]=mink(distanceSquared(:),min(8,numel(distanceSquared)));
for candidate=1:numel(values)
    if values(candidate)>4e-4; break; end
    [indexA,indexB]=ind2sub(size(distanceSquared),linearIndices(candidate));
    stateA=first.States{indexA}; stateB=second.States{indexB};
    relativeDistance=norm(stateA-stateB)/max([1,norm(stateA),norm(stateB)]);
    if relativeDistance<2e-3 && ...
            abs(a.freezeWeight(indexA)-b.freezeWeight(indexB))<2e-3
        intersects=true;
        return
    end
end
end
