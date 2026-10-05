function result = l6ns_atlas_state_sector_spectrum(setup,pathwayName,state,w,outputFile,modeCount)
% Resolve an arbitrary-state Jacobian spectrum by verified translation sectors.

if nargin<6 || isempty(modeCount); modeCount=16; end
J=l6ns_state_jacobian(setup,pathwayName,state,w);
rows={};
sectorModes=cell(8,1);
sectorIndex=0;
for rowWave=0:3
    for columnWave=0:1
        sectorIndex=sectorIndex+1;
        Q=l6ns_atlas_sector_basis(setup.Context.MapSize,rowWave,columnWave);
        reduced=Q'*J*Q;
        modes=l6ns_eigenpairs(reduced,min(modeCount,size(reduced,1)-2), ...
            'largestreal',setup.Config);
        [~,nearestIndex]=min(abs(modes.Lambda-1));
        sectorModes{sectorIndex}=struct('RowWave',rowWave, ...
            'ColumnWave',columnWave,'Lambda',modes.Lambda, ...
            'RightReduced',modes.Right,'LeftReduced',modes.Left);
        rows(end+1,:)={string(pathwayName),w,rowWave,columnWave, ... %#ok<AGROW>
            max(real(modes.Lambda)),sum(real(modes.Lambda)>1+1e-8), ...
            modes.Lambda(nearestIndex),abs(modes.Lambda(nearestIndex)-1), ...
            modes.ConditionNumber(nearestIndex)};
    end
end
summary=cell2table(rows,'VariableNames',{'pathway','freezeWeight', ...
    'rowWaveNumber','columnWaveNumber','maximumRealEigenvalue', ...
    'unstableDimensionInComputedSector','eigenvalueNearestOne', ...
    'distanceToOne','nearestModeConditionNumber'});
result=struct('Summary',summary,'SectorModes',{sectorModes},'State',state(:));
if nargin>=5 && ~isempty(outputFile)
    [outputDir,stem]=fileparts(outputFile);
    if ~exist(outputDir,'dir'); mkdir(outputDir); end
    save(outputFile,'result','-v7.3');
    writetable(summary,fullfile(outputDir,[stem '.tsv']), ...
        'FileType','text','Delimiter','\t');
end
end
