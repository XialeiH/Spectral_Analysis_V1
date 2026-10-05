function result = l6ns_atlas_symmetry_sectors(setup,outputDir,modeCount)
% Decompose the persistent Jacobians by row-10 and column-20 translations.

if nargin<3 || isempty(modeCount); modeCount=24; end
if ~exist(outputDir,'dir'); mkdir(outputDir); end
mapSize=setup.Context.MapSize;
if ~isequal(mapSize,[40 40])
    error('The current sector basis expects a 40-by-40 map.');
end
rows={};
sectors=cell(0,1);
pathways={'L6','Inhibition'};
for pathwayIndex=1:numel(pathways)
    pathwayName=pathways{pathwayIndex};
    critical=setup.(pathwayName);
    J=l6ns_state_jacobian(setup,pathwayName,setup.Context.FixedPoint,critical.W);
    for rowWave=0:3
        for columnWave=0:1
            Q=l6ns_atlas_sector_basis(mapSize,rowWave,columnWave);
            reduced=Q'*J*Q;
            modes=l6ns_eigenpairs(reduced,min(modeCount,size(reduced,1)-2), ...
                'largestreal',setup.Config);
            [~,criticalIndex]=min(abs(modes.Lambda-1));
            sectorIndex=numel(sectors)+1;
            sectors{sectorIndex}=struct('Pathway',string(pathwayName), ...
                'RowWave',rowWave,'ColumnWave',columnWave,'Basis',Q, ...
                'Lambda',modes.Lambda,'RightReduced',modes.Right, ...
                'LeftReduced',modes.Left);
            rows(end+1,:)={string(pathwayName),rowWave,columnWave, ... %#ok<AGROW>
                max(real(modes.Lambda)),sum(real(modes.Lambda)>1+1e-8), ...
                modes.Lambda(criticalIndex),min(abs(modes.Lambda-1)), ...
                modes.ConditionNumber(criticalIndex)};
        end
    end
end
summary=cell2table(rows,'VariableNames',{'pathway','rowWaveNumber', ...
    'columnWaveNumber','maximumRealEigenvalue','unstableDimensionInComputedSector', ...
    'eigenvalueNearestOne','distanceToOne','nearestModeConditionNumber'});
result=struct('Summary',summary,'Sectors',{sectors}, ...
    'RowTranslation',10,'ColumnTranslation',20);
writetable(summary,fullfile(outputDir,'persistent_symmetry_sector_spectra.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputDir,'persistent_symmetry_sector_spectra.mat'),'result','-v7.3');
end
