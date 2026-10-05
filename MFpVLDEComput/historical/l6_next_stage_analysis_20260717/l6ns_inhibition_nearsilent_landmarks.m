function tableOut = l6ns_inhibition_nearsilent_landmarks(setup,branchResult,outputDir)
% Evaluate stability at the saddle-node sample and both w_I=1 branch crossings.

fixed = setup.Context.FixedPoint(:);
phi = @(x,w)l6ns_phi(x,0,setup.Context,[1 1],1-w);
tableBranch = branchResult.Table;
states = branchResult.States;
[~,foldIndex] = min(tableBranch.freezeWeight);
foldIndices = unique(max(1,min(height(tableBranch),foldIndex+(-1:1))));
foldIndices = foldIndices(:);
landmarkStates = states(foldIndices);
landmarkWeights = tableBranch.freezeWeight(foldIndices);
landmarkNames = repmat("fold_neighbor",numel(foldIndices),1);
landmarkNames(foldIndices==foldIndex) = "minimum_weight_sample";

targetW = 1;
for index = 1:height(tableBranch)-1
    firstW = tableBranch.freezeWeight(index);
    secondW = tableBranch.freezeWeight(index+1);
    if (firstW-targetW)*(secondW-targetW)>0 || firstW==secondW
        continue
    end
    fraction = (targetW-firstW)/(secondW-firstW);
    guess = states{index}+fraction*(states{index+1}-states{index});
    [state,diagnostic] = l6ns_fixed_point_polish(phi,fixed,targetW,guess);
    if diagnostic.converged
        landmarkStates{end+1,1} = state; %#ok<AGROW>
        landmarkWeights(end+1,1) = targetW; %#ok<AGROW>
        landmarkNames(end+1,1) = "w1_crossing"; %#ok<AGROW>
    end
end

rows = {};
for index = 1:numel(landmarkStates)
    state = landmarkStates{index};
    w = landmarkWeights(index);
    residualNorm = norm(phi(state,w)-state);
    modes = l6ns_numeric_leading_modes(@(x)phi(x,w),state,4,setup.Config);
    rows(end+1,:) = {index,landmarkNames(index),w,norm(state-fixed), ... %#ok<AGROW>
        norm(state-fixed)/norm(fixed),max(real(modes.Lambda)), ...
        min(state),max(state),residualNorm,residualNorm/max(1,norm(fixed))};
end
tableOut = cell2table(rows,'VariableNames',{'landmarkIndex','landmark', ...
    'freezeWeight','distanceFromPersistent','relativeDistance','maxRealLambda', ...
    'minimumState','maximumState','residualNorm','relativeResidual'});
if ~exist(outputDir,'dir'); mkdir(outputDir); end
writetable(tableOut,fullfile(outputDir,'inhibition_nearsilent_landmarks.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputDir,'inhibition_nearsilent_landmarks.mat'), ...
    'tableOut','landmarkStates','-v7.3');
end
