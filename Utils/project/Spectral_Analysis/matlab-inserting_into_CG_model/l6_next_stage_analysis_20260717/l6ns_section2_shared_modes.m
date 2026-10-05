function result = l6ns_section2_shared_modes(cfg, data, section1)
% Section 2: robust eigen, singular, Fourier, and transient subspaces.

sectionDir = fullfile(cfg.OutputRoot, 'section2_shared_modes');
if ~exist(sectionDir, 'dir'); mkdir(sectionDir); end

wGrid = unique([cfg.SubspaceWGrid(:); section1.Summary.directWc]);
criteria = {'largest_real','closest_to_one','largest_modulus'};
dimensions = cfg.SubspaceDimensions;
subspaces = struct();
pairRows = {};

for criterionIndex = 1:numel(criteria)
    criterion = criteria{criterionIndex};
    for dimensionIndex = 1:numel(dimensions)
        k = dimensions(dimensionIndex);
        bases = cell(numel(wGrid), 1);
        for wi = 1:numel(wGrid)
            jacobian = data.A + (1 - wGrid(wi)) * data.B;
            selector = local_selector(criterion);
            invariant = l6ns_invariant_subspace(jacobian,k,selector,cfg);
            bases{wi} = invariant.Basis;
        end
        key = sprintf('%s_k%d', criterion, k);
        subspaces.(key) = bases;

        medianCosine = nan(numel(wGrid));
        minimumCosine = nan(numel(wGrid));
        grassmann = nan(numel(wGrid));
        for i = 1:numel(wGrid)
            for j = 1:numel(wGrid)
                cosine = svd(bases{i}' * bases{j});
                medianCosine(i,j) = median(cosine);
                minimumCosine(i,j) = min(cosine);
                grassmann(i,j) = sqrt(max(0, k - norm(bases{i}' * bases{j}, 'fro')^2));
                pairRows(end+1,:) = {criterion, k, wGrid(i), wGrid(j), ...
                    medianCosine(i,j), minimumCosine(i,j), grassmann(i,j)}; %#ok<AGROW>
            end
        end
        local_plot_subspace_heatmaps(wGrid, medianCosine, minimumCosine, grassmann, ...
            criterion, k, sectionDir);
    end
end

pairTable = cell2table(pairRows, 'VariableNames', ...
    {'criterion','dimension','w1','w2','medianCosine','minimumCosine','grassmannDistance'});
writetable(pairTable, fullfile(sectionDir, 'pairwise_subspace_geometry.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

% Consensus leading-real subspaces and basis-invariant leverage/Fourier maps.
consensusRows = {};
consensus = struct();
fourierTables = cell(numel(dimensions),1);
for dimensionIndex = 1:numel(dimensions)
    k = dimensions(dimensionIndex);
    bases = subspaces.(sprintf('largest_real_k%d', k));
    qStack = cell2mat(bases') / sqrt(numel(bases));
    [u,s,~] = svd(qStack, 'econ');
    qCommon = u(:,1:k);
    persistence = diag(s).^2;
    leakage = nan(numel(wGrid), 1);
    for wi = 1:numel(wGrid)
        jacobian = data.A + (1 - wGrid(wi)) * data.B;
        action = jacobian * qCommon;
        leakage(wi) = norm(action - qCommon * (qCommon' * action), 'fro') / ...
            max(norm(action, 'fro'), eps);
        consensusRows(end+1,:) = {k, wGrid(wi), leakage(wi), ...
            min(persistence(1:k)), median(persistence(1:k)), max(persistence(1:k))}; %#ok<AGROW>
    end
    consensus.(sprintf('k%d',k)) = struct('Basis',qCommon,'Persistence',persistence,'Leakage',leakage);
    fourierTables{dimensionIndex} = local_fourier_metrics(qCommon,data,k);
    local_plot_consensus_maps(qCommon, persistence, data, k, sectionDir);
end
consensusTable = cell2table(consensusRows, 'VariableNames', ...
    {'dimension','w','leakage','minimumPersistence','medianPersistence','maximumPersistence'});
writetable(consensusTable, fullfile(sectionDir, 'consensus_subspace.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
fourierTable = vertcat(fourierTables{:});
writetable(fourierTable,fullfile(sectionDir,'consensus_fourier_metrics.tsv'), ...
    'FileType','text','Delimiter','\t');

% Singular input-output channels and overlap with critical eigenmodes.
singularW = unique([0; 1; section1.Summary.directWc]);
singularRows = {};
singularResults = struct();
for wi = 1:numel(singularW)
    w = singularW(wi);
    jacobian = data.A + (1 - w) * data.B;
    topCount = min(cfg.SingularModeCount, size(jacobian,1)-2);
    [u,s,v] = svds(jacobian, topCount, 'largest');
    sigma = diag(s);
    nullAvailable = true;
    try
        nullCount = min(8,size(jacobian,1)-2);
        [uNull,sNull,vNull] = svds(jacobian,nullCount,'smallest');
        nullSigma = diag(sNull);
        [nullSigma,nullOrder] = sort(nullSigma,'ascend');
        uNull = uNull(:,nullOrder);
        vNull = vNull(:,nullOrder);
    catch exception
        warning('L6NS:NearNullJacobian', ...
            'Near-null singular projector unavailable at w=%.6g: %s',w,exception.message);
        uNull = [];
        vNull = [];
        nullSigma = [];
        nullAvailable = false;
    end
    fullSpectrumComputed = cfg.ComputeFullSingularSpectrum;
    if fullSpectrumComputed
        allSigma = svd(full(jacobian));
        allSigma = sort(allSigma,'descend');
        totalEnergy = sum(allSigma.^2);
        numericalRanks = arrayfun(@(threshold) ...
            sum(allSigma > allSigma(1)*threshold),[1e-6 1e-8 1e-10 1e-12]);
        cumulativeEnergy = cumsum(allSigma.^2) / max(totalEnergy,eps);
        spectrumTable = table((1:numel(allSigma))',allSigma,allSigma/allSigma(1), ...
            cumulativeEnergy,'VariableNames',{'index','singularValue', ...
            'relativeSingularValue','cumulativeEnergyFraction'});
        writetable(spectrumTable,fullfile(sectionDir, ...
            sprintf('full_singular_spectrum_%s.tsv',local_w_field(w))), ...
            'FileType','text','Delimiter','\t');
    else
        allSigma = sigma;
        totalEnergy = norm(jacobian,'fro')^2;
        numericalRanks = nan(1,4);
        cumulativeEnergy = cumsum(allSigma.^2) / max(totalEnergy,eps);
    end
    stableRank = totalEnergy / allSigma(1)^2;
    consecutiveGap = allSigma(1:end-1) ./ ...
        max(allSigma(2:end),eps*allSigma(1));
    [largestGapRatio,largestGapIndex] = max(consecutiveGap);
    smallestRelativeSingularValue = allSigma(end)/allSigma(1);
    criticalModes = l6ns_eigenpairs(jacobian, min(20,topCount), 'largestreal', cfg);
    qEig = orth(criticalModes.Right);
    overlapRight = svd(qEig' * v(:,1:min(20,size(v,2))));
    overlapLeft = svd(qEig' * u(:,1:min(20,size(u,2))));
    singularRows(end+1,:) = {w,fullSpectrumComputed,nullAvailable,sigma(1), ...
        local_first_or_nan(nullSigma),stableRank, ...
        numericalRanks(1),numericalRanks(2),numericalRanks(3),numericalRanks(4), ...
        sum(sigma > sigma(1)*1e-6), sum(sigma > sigma(1)*1e-8), ...
        largestGapRatio,largestGapIndex,smallestRelativeSingularValue, ...
        cumulativeEnergy(min(20,numel(cumulativeEnergy))), ...
        median(overlapRight), min(overlapRight), median(overlapLeft), min(overlapLeft)}; %#ok<AGROW>
    singularResults.(local_w_field(w)) = struct('U',u,'Sigma',sigma,'V',v, ...
        'NullU',uNull,'NullSigma',nullSigma,'NullV',vNull, ...
        'AllSigma',allSigma,'CumulativeEnergy',cumulativeEnergy, ...
        'FullSpectrumComputed',fullSpectrumComputed);
    local_plot_singular_channels(u(:,1),v(:,1),allSigma,cumulativeEnergy,data,w,sectionDir);
    if nullAvailable
        local_plot_near_null_projector(uNull,vNull,data,w,sectionDir);
    end
end
singularTable = cell2table(singularRows, 'VariableNames', ...
    {'w','fullSpectrumComputed','nearNullProjectorAvailable', ...
    'largestSingularValue','smallestComputedSingularValue','stableRank', ...
    'numericalRankRelative1e6','numericalRankRelative1e8', ...
    'numericalRankRelative1e10','numericalRankRelative1e12', ...
    'computedTopCountRelative1e6','computedTopCountRelative1e8', ...
    'largestConsecutiveGapRatio','largestGapIndex','smallestRelativeSingularValue', ...
    'computedTop20EnergyFraction','medianEigRightSingularCosine','minimumEigRightSingularCosine', ...
    'medianEigLeftSingularCosine','minimumEigLeftSingularCosine'});
writetable(singularTable, fullfile(sectionDir, 'singular_subspace_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

% Continuous-time transient modes for the original, critical, and fixed-L6 systems.
transientResults = struct();
transientRows = {};
n=data.PopulationSize;
fixed=data.FixedPoint;
populationScale=[norm(fixed(1:n))/sqrt(n),norm(fixed(n+(1:n)))/sqrt(n), ...
    norm(fixed(2*n+(1:n)))/sqrt(n)];
populationScale=max(populationScale,1e-8);
transientWeights=[ones(n,1)/populationScale(1);ones(n,1)/populationScale(2); ...
    ones(n,1)/populationScale(3)];
for wi = 1:numel(singularW)
    w = singularW(wi);
    jacobian = data.A + (1 - w) * data.B;
    transient = l6ns_continuous_transient(jacobian,cfg.TransientTimeGrid,cfg,transientWeights);
    criticalModes = l6ns_eigenpairs(jacobian, 12, 'largestreal', cfg);
    qEig = orth(transientWeights.*criticalModes.Right);
    weightedInput=transientWeights.*transient.OptimalInput;
    weightedOutput=transientWeights.*transient.PeakOutput;
    overlapInput = norm(qEig'*weightedInput)^2/max(norm(weightedInput)^2,eps);
    overlapOutput = norm(qEig'*weightedOutput)^2/max(norm(weightedOutput)^2,eps);
    transientRows(end+1,:) = {w,transient.NormType,populationScale(1),populationScale(2), ...
        populationScale(3),transient.PeakTime,transient.PeakGain,overlapInput,overlapOutput}; %#ok<AGROW>
    transientResults.(local_w_field(w)) = transient;
    local_plot_transient(transient, data, w, sectionDir);
end
transientTable = cell2table(transientRows, 'VariableNames', ...
    {'w','normType','S_scale','C_scale','I_scale','peakTime','peakGain', ...
    'criticalOverlapInput','criticalOverlapOutput'});
writetable(transientTable, fullfile(sectionDir, 'continuous_transient_summary.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

result = struct();
result.Pairwise = pairTable;
result.Consensus = consensusTable;
result.ConsensusObjects = consensus;
result.Fourier = fourierTable;
result.Singular = singularTable;
result.SingularObjects = singularResults;
result.Transient = transientTable;
result.TransientObjects = transientResults;
save(fullfile(sectionDir, 'section2_result.mat'), 'result', '-v7.3');
end

function tableOut = local_fourier_metrics(q,data,k)
n=data.PopulationSize;
mapS=local_fourier_stack(q(1:n,:),data.MapSize);
mapC=local_fourier_stack(q(n+(1:n),:),data.MapSize);
mapI=local_fourier_stack(q(2*n+(1:n),:),data.MapSize);
mapE=(1-data.ExcitatoryCWeight)*reshape(q(1:n,:),[data.MapSize size(q,2)])+ ...
    data.ExcitatoryCWeight*reshape(q(n+(1:n),:),[data.MapSize size(q,2)]);
mapE=local_fourier_stack(reshape(mapE,[],size(q,2)),data.MapSize);
stacks={mapS,mapC,mapI,mapE};

[ny,nx,~]=size(mapS);
cy=floor(ny/2)+1; cx=floor(nx/2)+1;
[gridY,gridX]=ndgrid(((1:ny)-cy)/ny,((1:nx)-cx)/nx);
radial=hypot(gridX,gridY);
orientation=mod(atan2d(gridY,gridX),180);

ePower=sum(abs(mapE).^2,3); ePower(cy,cx)=0;
[~,dominantIndex]=max(ePower(:));
[dominantY,dominantX]=ind2sub([ny nx],dominantIndex);
phaseSC=angle(sum(conj(mapS(dominantY,dominantX,:)).*mapC(dominantY,dominantX,:),'all'))*180/pi;
phaseSI=angle(sum(conj(mapS(dominantY,dominantX,:)).*mapI(dominantY,dominantX,:),'all'))*180/pi;
phaseCI=angle(sum(conj(mapC(dominantY,dominantX,:)).*mapI(dominantY,dominantX,:),'all'))*180/pi;

rows=cell(4,1);
labels={'S','C','I','E'};
for population=1:4
    power=sum(abs(stacks{population}).^2,3);
    power(cy,cx)=0;
    total=sum(power,'all');
    [~,index]=max(power(:));
    meanFrequency=sum(radial.*power,'all')/max(total,eps);
    bandwidth=sqrt(sum((radial-meanFrequency).^2.*power,'all')/max(total,eps));
    kernelOverlap=NaN;
    if ~isempty(data.L6Kernel)
        kernelPower=abs(fftshift(fft2(data.L6Kernel,ny,nx))).^2;
        kernelPower(cy,cx)=0;
        kernelOverlap=sum(power.*kernelPower,'all')/ ...
            max(norm(power(:))*norm(kernelPower(:)),eps);
    end
    rows{population}={k,labels{population},radial(index),orientation(index), ...
        bandwidth,kernelOverlap,phaseSC,phaseSI,phaseCI};
end
tableOut=cell2table(vertcat(rows{:}),'VariableNames', ...
    {'dimension','population','dominantSpatialFrequency','dominantOrientationDeg', ...
    'frequencyBandwidth','l6KernelPowerOverlap','phaseSC_deg','phaseSI_deg','phaseCI_deg'});
end

function fourierStack = local_fourier_stack(vectors,mapSize)
fourierStack=zeros([mapSize size(vectors,2)]);
for column=1:size(vectors,2)
    fourierStack(:,:,column)=fftshift(fft2(reshape(vectors(:,column),mapSize)));
end
end

function selector = local_selector(criterion)
switch criterion
    case 'largest_real'
        selector = 'largestreal';
    case 'closest_to_one'
        selector = 1;
    case 'largest_modulus'
        selector = 'largestabs';
    otherwise
        error('Unknown subspace criterion %s.', criterion);
end
end

function local_plot_subspace_heatmaps(wGrid, medianCosine, minimumCosine, grassmann, criterion, k, outputDir)
fig = figure('Visible','off','Color','w','Position',[80 80 1500 430]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
values = {medianCosine, minimumCosine, grassmann};
titles = {'median cosine','minimum cosine','Grassmann distance'};
for i = 1:3
    nexttile; imagesc(wGrid,wGrid,values{i}); axis xy image; colorbar;
    xlabel('w_j'); ylabel('w_i'); title(titles{i});
end
sgtitle(sprintf('%s subspace, k=%d', strrep(criterion,'_',' '), k));
l6ns_save_figure(fig, outputDir, sprintf('subspace_%s_k%d',criterion,k));
close(fig);
end

function local_plot_consensus_maps(q, persistence, data, k, outputDir)
n = data.PopulationSize;
leverageS = sum(abs(q(1:n,:)).^2,2);
leverageC = sum(abs(q(n+(1:n),:)).^2,2);
leverageI = sum(abs(q(2*n+(1:n),:)).^2,2);
qE=(1-data.ExcitatoryCWeight)*q(1:n,:)+data.ExcitatoryCWeight*q(n+(1:n),:);
leverageE = sum(abs(qE).^2,2);
maps = {leverageS,leverageC,leverageI,leverageE};
populationLabels={'S','C','I','E'};
fig = figure('Visible','off','Color','w','Position',[80 80 1450 760]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for i = 1:4
    nexttile; imagesc(reshape(maps{i},data.MapSize)); axis image off; colorbar;
    title(sprintf('%s leverage',populationLabels{i}));
end
for population = 1:4
    nexttile;
    power = local_subspace_fourier_power(q, population, data);
    imagesc(fftshift(power)); axis image off; colorbar;
    title(sprintf('%s Fourier power',populationLabels{population}));
end
sgtitle(sprintf('Consensus critical subspace k=%d, persistence %.3f--%.3f', ...
    k, min(persistence(1:k)), max(persistence(1:k))));
l6ns_save_figure(fig, outputDir, sprintf('consensus_leverage_fourier_k%d',k));
close(fig);
end

function power = local_subspace_fourier_power(q, population, data)
n = data.PopulationSize;
power = zeros(data.MapSize);
for column = 1:size(q,2)
    switch population
        case 1
            map = reshape(q(1:n,column),data.MapSize);
        case 2
            map = reshape(q(n+(1:n),column),data.MapSize);
        case 3
            map = reshape(q(2*n+(1:n),column),data.MapSize);
        case 4
            map = (1-data.ExcitatoryCWeight)*reshape(q(1:n,column),data.MapSize) + ...
                data.ExcitatoryCWeight*reshape(q(n+(1:n),column),data.MapSize);
    end
    power = power + abs(fft2(map)).^2;
end
end

function local_plot_singular_channels(u, v, sigma, cumulativeEnergy, data, w, outputDir)
n = data.PopulationSize;
populationLabels={'S','C','I','E'};
fig = figure('Visible','off','Color','w','Position',[80 80 1500 850]);
tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
for row = 1:2
    if row==1
        vector=v; label='input v_1';
    else
        vector=u; label='output u_1';
    end
    maps = {reshape(real(vector(1:n)),data.MapSize), ...
        reshape(real(vector(n+(1:n))),data.MapSize), ...
        reshape(real(vector(2*n+(1:n))),data.MapSize)};
    maps{4} = (1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit = max(cellfun(@(x) max(abs(x(:))),maps));
    for col = 1:4
        nexttile; imagesc(maps{col}); axis image off; colorbar;
        if limit>0; clim([-limit limit]); end
        title(sprintf('%s %s',label,populationLabels{col}));
    end
end
nexttile([1 2]); semilogy(sigma,'-o'); grid on; xlabel('index'); ylabel('\sigma'); title('leading singular values');
nexttile([1 2]); plot(cumulativeEnergy,'-o'); grid on; xlabel('index'); ylabel('fraction of ||J||_F^2'); title('cumulative energy');
sgtitle(sprintf('Singular input-output channels, w=%.5g',w));
l6ns_save_figure(fig, outputDir, sprintf('singular_channels_%s',local_w_field(w)));
close(fig);
end

function local_plot_near_null_projector(uNull,vNull,data,w,outputDir)
n=data.PopulationSize;
populationLabels={'S','C','I','E'};
fig=figure('Visible','off','Color','w','Position',[80 80 1450 700]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for row=1:2
    if row==1; basis=vNull; side='input'; else; basis=uNull; side='output'; end
    maps={sum(abs(basis(1:n,:)).^2,2),sum(abs(basis(n+(1:n),:)).^2,2), ...
        sum(abs(basis(2*n+(1:n),:)).^2,2)};
    excitatory=(1-data.ExcitatoryCWeight)*basis(1:n,:)+ ...
        data.ExcitatoryCWeight*basis(n+(1:n),:);
    maps{4}=sum(abs(excitatory).^2,2);
    for column=1:4
        nexttile; imagesc(reshape(maps{column},data.MapSize)); axis image off; colorbar;
        title(sprintf('%s %s',side,populationLabels{column}));
    end
end
sgtitle(sprintf('Near-null singular projector leverage, w=%.5g',w));
l6ns_save_figure(fig,outputDir,sprintf('near_null_projector_%s',local_w_field(w)));
close(fig);
end

function value=local_first_or_nan(values)
if isempty(values); value=NaN; else; value=values(1); end
end

function local_plot_transient(transient, data, w, outputDir)
n = data.PopulationSize;
populationLabels={'S','C','I','E'};
fig = figure('Visible','off','Color','w','Position',[80 80 1500 760]);
tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
nexttile([1 1]); plot(transient.TimeGrid,transient.Gain,'-o'); grid on; xlabel('time'); ylabel('||e^{t(J-I)}||_2'); title('gain');
for row = 1:2
    if row==1
        vector=transient.OptimalInput; label='optimal input';
    else
        vector=transient.PeakOutput; label='peak output';
    end
    maps = {reshape(real(vector(1:n)),data.MapSize), ...
        reshape(real(vector(n+(1:n))),data.MapSize), ...
        reshape(real(vector(2*n+(1:n))),data.MapSize)};
    maps{4} = (1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    startCol = 2;
    for col = startCol:5
        tileIndex = (row-1)*5 + col;
        nexttile(tileIndex); imagesc(maps{col-1}); axis image off; colorbar;
        title(sprintf('%s %s',label,populationLabels{col-1}));
    end
end
sgtitle(sprintf('Continuous-time transient, w=%.5g, peak t=%.3g, gain=%.4g', ...
    w,transient.PeakTime,transient.PeakGain));
l6ns_save_figure(fig, outputDir, sprintf('continuous_transient_%s',local_w_field(w)));
close(fig);
end

function field = local_w_field(w)
field = sprintf('w_%+.6f',w);
field = strrep(field,'+','p');
field = strrep(field,'-','m');
field = strrep(field,'.','p');
end
