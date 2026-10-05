function inventory = l6ns_atlas_symmetry_inventory(setup,samples,outputDir,tolerance)
% Audit translations/reflections and label each sampled state's isotropy.

if nargin<4 || isempty(tolerance); tolerance=1e-8; end
if ~exist(outputDir,'dir'); mkdir(outputDir); end
mapSize=setup.Context.MapSize;
rowShifts=0:mapSize(1)/4:mapSize(1)-mapSize(1)/4;
columnShifts=0:mapSize(2)/4:mapSize(2)-mapSize(2)/4;
rows={};
candidateIndex=0;
for reflectRows=[false true]
    for reflectColumns=[false true]
        for rowShift=rowShifts
            for columnShift=columnShifts
                candidateIndex=candidateIndex+1;
                errors=nan(numel(samples),1);
                invariance=nan(numel(samples),1);
                for sampleIndex=1:numel(samples)
                    sample=samples{sampleIndex};
                    phi=local_phi(setup,sample.Pathway,sample.W);
                    state=sample.State(:);
                    transformed=l6ns_atlas_transform_state(state,mapSize,rowShift, ...
                        columnShift,reflectRows,reflectColumns);
                    residual=phi(state)-state;
                    transformedResidual=phi(transformed)-transformed;
                    expected=l6ns_atlas_transform_state(residual,mapSize,rowShift, ...
                        columnShift,reflectRows,reflectColumns);
                    errors(sampleIndex)=norm(transformedResidual-expected)/(1+norm(residual));
                    invariance(sampleIndex)=norm(transformed-state)/max(1,norm(state));
                end
                label=sprintf('r%d_c%d_fr%d_fc%d',rowShift,columnShift, ...
                    reflectRows,reflectColumns);
                rows(end+1,:)={candidateIndex,string(label),rowShift,columnShift, ... %#ok<AGROW>
                    reflectRows,reflectColumns,max(errors),median(errors), ...
                    max(invariance),max(errors)<tolerance};
            end
        end
    end
end
candidateTable=cell2table(rows,'VariableNames',{'candidateIndex','label', ...
    'rowShift','columnShift','reflectRows','reflectColumns', ...
    'maximumEquivarianceError','medianEquivarianceError', ...
    'maximumSampleInvarianceError','verifiedEquationSymmetry'});
verified=candidateTable(candidateTable.verifiedEquationSymmetry,:);

isotropyRows={};
for sampleIndex=1:numel(samples)
    sample=samples{sampleIndex};
    state=sample.State(:);
    invariantLabels=strings(0,1);
    for symmetryIndex=1:height(verified)
        transformed=l6ns_atlas_transform_state(state,mapSize, ...
            verified.rowShift(symmetryIndex),verified.columnShift(symmetryIndex), ...
            verified.reflectRows(symmetryIndex),verified.reflectColumns(symmetryIndex));
        errorValue=norm(transformed-state)/max(1,norm(state));
        if errorValue<tolerance
            invariantLabels(end+1,1)=verified.label(symmetryIndex); %#ok<AGROW>
        end
    end
    isotropySize=numel(invariantLabels);
    orbitSize=height(verified)/max(1,isotropySize);
    isotropyRows(end+1,:)={sampleIndex,string(sample.Label),string(sample.Pathway), ... %#ok<AGROW>
        sample.W,isotropySize,orbitSize,strjoin(invariantLabels,',')};
end
isotropyTable=cell2table(isotropyRows,'VariableNames',{'sampleIndex','sampleLabel', ...
    'pathway','freezeWeight','isotropySize','orbitSizeWithinVerifiedGroup','isotropyLabels'});

writetable(candidateTable,fullfile(outputDir,'symmetry_candidate_audit.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(isotropyTable,fullfile(outputDir,'sample_isotropy.tsv'), ...
    'FileType','text','Delimiter','\t');
inventory=struct('Candidates',candidateTable,'VerifiedGroup',verified, ...
    'Isotropy',isotropyTable,'Tolerance',tolerance);
save(fullfile(outputDir,'symmetry_inventory.mat'),'inventory','-v7.3');
end

function phi=local_phi(setup,pathwayName,w)
if strcmpi(pathwayName,'L6')
    phi=@(x)l6ns_phi(x,w,setup.Context,[1 1],1);
else
    phi=@(x)l6ns_phi(x,0,setup.Context,[1 1],1-w);
end
end
