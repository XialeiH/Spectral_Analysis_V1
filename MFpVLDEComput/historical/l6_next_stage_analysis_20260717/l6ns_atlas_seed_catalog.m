function catalog = l6ns_atlas_seed_catalog(setup,aggregateFile,nearSilentFile,outputDir)
% Collect every known census root, branch intersection, fold, and stable seed.

if ~exist(outputDir,'dir'); mkdir(outputDir); end
loaded=load(aggregateFile,'result');
aggregate=loaded.result;
states={};
rows={};

for itemIndex=1:numel(aggregate.Multistart)
    item=aggregate.Multistart{itemIndex};
    for rootIndex=1:numel(item.RootStates)
        [states,rows]=local_add(states,rows,setup,item.Name,item.FreezeWeight, ...
            item.RootStates{rootIndex},sprintf('multistart_%d_root_%d',itemIndex,rootIndex), ...
            'census_root');
    end
end

pathways={'L6','Inhibition'};
for pathwayIndex=1:numel(pathways)
    pathwayName=pathways{pathwayIndex};
    branches=aggregate.Branch.(pathwayName);
    targetWeights=unique(aggregate.RootCensus.freezeWeight( ...
        aggregate.RootCensus.pathway==string(pathwayName)));
    for branchIndex=1:numel(branches)
        branch=branches{branchIndex};
        for targetIndex=1:numel(targetWeights)
            targetW=targetWeights(targetIndex);
            crossings=local_branch_crossings(setup,pathwayName,branch,targetW);
            for crossingIndex=1:numel(crossings)
                [states,rows]=local_add(states,rows,setup,pathwayName,targetW, ...
                    crossings{crossingIndex},sprintf('branch_%d_intersection_%d', ...
                    branchIndex,crossingIndex),'branch_intersection');
            end
        end
        extrema=local_branch_extrema(branch);
        for extremaIndex=1:numel(extrema)
            index=extrema(extremaIndex);
            [states,rows]=local_add(states,rows,setup,pathwayName, ...
                branch.Table.freezeWeight(index),branch.States{index}, ...
                sprintf('branch_%d_fold_seed_%d',branchIndex,index),'fold_seed');
        end
    end
end

if exist(nearSilentFile,'file')
    loaded=load(nearSilentFile,'result');
    branch=loaded.result;
    crossings=local_branch_crossings(setup,'Inhibition',branch,1);
    for crossingIndex=1:numel(crossings)
        [states,rows]=local_add(states,rows,setup,'Inhibition',1, ...
            crossings{crossingIndex},sprintf('nearsilent_w1_crossing_%d',crossingIndex), ...
            'nearsilent_family');
    end
    extrema=local_branch_extrema(branch);
    for extremaIndex=1:numel(extrema)
        index=extrema(extremaIndex);
        [states,rows]=local_add(states,rows,setup,'Inhibition', ...
            branch.Table.freezeWeight(index),branch.States{index}, ...
            sprintf('nearsilent_fold_seed_%d',index),'nearsilent_fold_seed');
    end
end

tableOut=cell2table(rows,'VariableNames',{'seedId','pathway','freezeWeight', ...
    'label','seedClass','distanceFromPersistent','minimumState','maximumState'});
catalog=struct('Table',tableOut,'States',{states});
writetable(tableOut,fullfile(outputDir,'atlas_seed_catalog.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputDir,'atlas_seed_catalog.mat'),'catalog','-v7.3');
end

function crossings=local_branch_crossings(setup,pathwayName,branch,targetW)
w=branch.Table.freezeWeight;
crossings={};
for index=1:numel(w)-1
    if (w(index)-targetW)*(w(index+1)-targetW)>0
        continue
    end
    denominator=w(index+1)-w(index);
    if abs(denominator)<eps
        fraction=0.5;
    else
        fraction=(targetW-w(index))/denominator;
    end
    if fraction<-1e-10 || fraction>1+1e-10; continue; end
    guess=(1-fraction)*branch.States{index}+fraction*branch.States{index+1};
    phi=local_phi(setup,pathwayName);
    [state,diagnostic]=l6ns_fixed_point_newton(phi,setup.Context.FixedPoint,targetW,guess);
    if diagnostic.converged
        crossings{end+1,1}=state; %#ok<AGROW>
    end
end
end

function indices=local_branch_extrema(branch)
w=branch.Table.freezeWeight(:);
smoothed=movmean(w,5,'Endpoints','shrink');
direction=sign(diff(smoothed));
direction(direction==0)=1;
candidates=find(direction(1:end-1).*direction(2:end)<0)+1;
window=4;
minimumProminence=5e-4;
keep=false(size(candidates));
for candidateIndex=1:numel(candidates)
    index=candidates(candidateIndex);
    if index<=window || index>numel(w)-window; continue; end
    prominence=min(abs(smoothed(index)-smoothed(index-window)), ...
        abs(smoothed(index)-smoothed(index+window)));
    keep(candidateIndex)=prominence>=minimumProminence;
end
indices=candidates(keep);
end

function [states,rows]=local_add(states,rows,setup,pathwayName,w,state,label,seedClass)
state=state(:);
for index=1:numel(states)
    existing=rows(index,:);
    if existing{2}==string(pathwayName) && abs(existing{3}-w)<1e-9 && ...
            norm(states{index}-state)/max(1,norm(state))<1e-7
        return
    end
end
seedId=numel(states)+1;
states{seedId,1}=state;
fixed=setup.Context.FixedPoint(:);
rows(seedId,:)={seedId,string(pathwayName),w,string(label),string(seedClass), ...
    norm(state-fixed),min(state),max(state)};
end

function phi=local_phi(setup,pathwayName)
if strcmpi(pathwayName,'L6')
    phi=@(x,w)l6ns_phi(x,w,setup.Context,[1 1],1);
else
    phi=@(x,w)l6ns_phi(x,0,setup.Context,[1 1],1-w);
end
end
