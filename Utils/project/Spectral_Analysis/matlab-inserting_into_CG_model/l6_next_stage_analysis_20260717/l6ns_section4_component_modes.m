function result = l6ns_section4_component_modes(cfg, data, section1)
% Section 4: component-mode atlas and adjoint L6 attribution.

sectionDir = fullfile(cfg.OutputRoot,'section4_component_modes');
if ~exist(sectionDir,'dir'); mkdir(sectionDir); end

operators = {'A_base','B_L6'};
matrices = {data.A,data.B};
categories = {'positive_edge','negative_edge','near_zero'};
% B has an exact nullspace, so a literal zero shift makes eigs' sparse
% shift-and-invert factorization singular. The tiny offset still selects
% the near-zero sector without changing the interpreted scale.
selectors = {'largestreal','smallestreal',1e-8};
atlasRows = {};
modeObjects = struct();

criticalR = section1.CriticalRight;
criticalL = section1.CriticalLeft;
for operatorIndex = 1:numel(operators)
    operatorName = operators{operatorIndex};
    matrix = matrices{operatorIndex};
    for categoryIndex = 1:numel(categories)
        category = categories{categoryIndex};
        modes = l6ns_eigenpairs(matrix,cfg.ComponentAtlasSize,selectors{categoryIndex},cfg);
        key = sprintf('%s_%s',operatorName,category);
        modeObjects.(key) = modes;
        for modeIndex = 1:numel(modes.Lambda)
            r = modes.Right(:,modeIndex);
            fingerprint = local_fingerprint(r,data);
            atlasRows(end+1,:) = {operatorName,category,modeIndex, ...
                real(modes.Lambda(modeIndex)),imag(modes.Lambda(modeIndex)), ...
                modes.ConditionNumber(modeIndex),modes.RightResidual(modeIndex), ...
                fingerprint.SFraction,fingerprint.CFraction,fingerprint.IFraction, ...
                fingerprint.EFraction,fingerprint.ParticipationRatio, ...
                fingerprint.DominantFrequency,fingerprint.DominantOrientation, ...
                abs(criticalL'*matrix*r),abs(criticalR'*r)}; %#ok<AGROW>
        end
local_plot_mode_atlas(modes.Right(:,1:min(4,size(modes.Right,2))), ...
            modes.Lambda(1:min(4,numel(modes.Lambda))),data,operatorName,category,sectionDir);
        complexIndex=find(abs(imag(modes.Lambda))>1e-8,1,'first');
        if ~isempty(complexIndex)
            local_plot_complex_mode(modes.Right(:,complexIndex),modes.Lambda(complexIndex), ...
                data,operatorName,category,sectionDir);
        end
    end
end

% Leading singular input-output pairs and basis-invariant near-null sectors.
[singularTable,singularSummary,singularObjects] = local_component_singular_analysis( ...
    operators,matrices,cfg,data,sectionDir);
writetable(singularTable,fullfile(sectionDir,'component_singular_atlas.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(singularSummary,fullfile(sectionDir,'component_singular_summary.tsv'), ...
    'FileType','text','Delimiter','\t');

atlasTable = cell2table(atlasRows,'VariableNames',{'operator','category','mode', ...
    'lambdaReal','lambdaImag','conditionNumber','rightResidual', ...
    'S_energyFraction','C_energyFraction','I_energyFraction','E_energyFraction', ...
    'participationRatio','dominantSpatialFrequency','dominantSpatialOrientationDeg', ...
    'criticalAdjointAction','criticalRightOverlap'});
writetable(atlasTable,fullfile(sectionDir,'component_mode_atlas.tsv'), ...
    'FileType','text','Delimiter','\t');

% Compare positive-edge subspaces of A, B, and the full critical Jacobian.
k = cfg.PositiveSubspaceSize;
aPositive = l6ns_invariant_subspace(data.A,k,'largestreal',cfg);
bPositive = l6ns_invariant_subspace(data.B,k,'largestreal',cfg);
qA = aPositive.Basis;
qB = bPositive.Basis;
jCritical = data.A+(1-section1.Summary.directWc)*data.B;
jInvariant = l6ns_invariant_subspace(jCritical,k,'largestreal',cfg);
qJ = jInvariant.Basis;
subspaceTable = local_subspace_comparison(qA,qB,qJ);
writetable(subspaceTable,fullfile(sectionDir,'component_positive_subspaces.tsv'), ...
    'FileType','text','Delimiter','\t');

% Stabilizing/destabilizing actions on every full-system mode in the full
% run, with a reduced set retained for smoke validation.
if cfg.ComputeFullEigenSpectrum
    fullModes = local_full_eigenpairs(jCritical);
else
    fullModes = l6ns_eigenpairs(jCritical,cfg.AnalysisModeCount,'largestreal',cfg);
end
attributionRows = cell(numel(fullModes.Lambda),1);
alphaCritical = 1-section1.Summary.directWc;
actionAAll = data.A*fullModes.Right;
actionBAll = alphaCritical*(data.B*fullModes.Right);
for modeIndex = 1:numel(fullModes.Lambda)
    l = fullModes.Left(:,modeIndex);
    cA = l'*actionAAll(:,modeIndex);
    cB = l'*actionBAll(:,modeIndex);
    lambda = fullModes.Lambda(modeIndex);
    radialA = real(conj(lambda)/max(abs(lambda),eps)*cA);
    radialB = real(conj(lambda)/max(abs(lambda),eps)*cB);
    status = "stable";
    if real(lambda)>=1; status="unstable";
    elseif real(lambda)>=0.9; status="near_critical";
    end
    attributionRows{modeIndex} = {modeIndex,real(lambda),imag(lambda), ...
        real(cA),imag(cA),real(cB),imag(cB),radialA,radialB, ...
        fullModes.ConditionNumber(modeIndex),status};
end
attributionTable = cell2table(vertcat(attributionRows{:}),'VariableNames', ...
    {'mode','lambdaReal','lambdaImag','baseReal','baseImag','l6Real','l6Imag', ...
    'baseRadial','l6Radial','conditionNumber','status'});
writetable(attributionTable,fullfile(sectionDir,'full_mode_component_attribution.tsv'), ...
    'FileType','text','Delimiter','\t');
local_plot_attribution_scatter(attributionTable,sectionDir);

% How one component redirects modes of the other.
redirectAB = local_redirection(data.B,qA,'B_on_A');
redirectBA = local_redirection(data.A,qB,'A_on_B');
redirectionTable = [redirectAB;redirectBA];
writetable(redirectionTable,fullfile(sectionDir,'component_mode_redirection.tsv'), ...
    'FileType','text','Delimiter','\t');
local_plot_redirection(redirectionTable,sectionDir);
local_plot_cross_actions(qA(:,1),qB(:,1),data.A,data.B,data,sectionDir);

% Exact operator-level adjoint attribution of Re(l^*Br), localized at both
% output targets and input sources.
br = data.B*criticalR;
adjointInput = data.B'*criticalL;
targetDensity = real(conj(criticalL).*br);
sourceDensity = real(conj(adjointInput).*criticalR);
targetTotal = sum(targetDensity);
sourceTotal = sum(sourceDensity);
scalarTotal = real(criticalL'*data.B*criticalR);
adjointSummary = table(targetTotal,sourceTotal,scalarTotal, ...
    abs(targetTotal-scalarTotal),abs(sourceTotal-scalarTotal), ...
    'VariableNames',{'targetDensitySum','sourceDensitySum','scalarSensitivity', ...
    'targetClosureError','sourceClosureError'});
writetable(adjointSummary,fullfile(sectionDir,'adjoint_attribution_closure.tsv'), ...
    'FileType','text','Delimiter','\t');
local_plot_adjoint_density(targetDensity,sourceDensity,data,sectionDir);

[l6InputOperator,l6OutputOperator,factorResidual] = local_factor_l6_operator(data.B,data.PopulationSize);
z=l6InputOperator*criticalR;
p=l6OutputOperator'*criticalL;
n=data.PopulationSize;
targetChannelDensity={ ...
    real(conj(l6OutputOperator(1:n,:)'*criticalL(1:n)).*z), ...
    real(conj(l6OutputOperator(n+(1:n),:)'*criticalL(n+(1:n))).*z), ...
    real(conj(l6OutputOperator(2*n+(1:n),:)'*criticalL(2*n+(1:n))).*z)};
zSourceS=l6InputOperator(:,1:n)*criticalR(1:n);
zSourceC=l6InputOperator(:,n+(1:n))*criticalR(n+(1:n));
sourceChannelDensity={real(conj(p).*zSourceS),real(conj(p).*zSourceC)};
factorDensity=real(conj(p).*z);
factorSummary=table(factorResidual,sum(factorDensity),scalarTotal, ...
    sum(targetChannelDensity{1}),sum(targetChannelDensity{2}),sum(targetChannelDensity{3}), ...
    sum(sourceChannelDensity{1}),sum(sourceChannelDensity{2}), ...
    abs(sum(factorDensity)-scalarTotal), ...
    'VariableNames',{'factorReconstructionResidual','factorDensitySum','scalarSensitivity', ...
    'targetSContribution','targetCContribution','targetIContribution', ...
    'sourceSContribution','sourceCContribution','factorClosureError'});
writetable(factorSummary,fullfile(sectionDir,'physical_l6_factor_attribution.tsv'), ...
    'FileType','text','Delimiter','\t');
local_plot_l6_factor_density(targetChannelDensity,sourceChannelDensity,factorDensity,data,sectionDir);

result = struct('Atlas',atlasTable,'PositiveSubspaces',subspaceTable, ...
    'FullModeAttribution',attributionTable,'Redirection',redirectionTable, ...
    'AdjointSummary',adjointSummary,'ModeObjects',modeObjects, ...
    'PhysicalFactorSummary',factorSummary, ...
    'SingularAtlas',singularTable,'SingularSummary',singularSummary, ...
    'SingularObjects',singularObjects);
save(fullfile(sectionDir,'section4_result.mat'),'result','-v7.3');
end

function modes = local_full_eigenpairs(matrix)
[right,lambda,left] = eig(full(matrix),'vector');
count = numel(lambda);
conditionNumber = nan(count,1);
for index = 1:count
    r = right(:,index);
    r = r/max(norm(r),eps);
    l = left(:,index);
    overlap = l'*r;
    l = l/conj(overlap);
    right(:,index) = r;
    left(:,index) = l;
    conditionNumber(index) = norm(l);
end
[~,order] = sort(real(lambda),'descend');
modes = struct('Lambda',lambda(order),'Right',right(:,order), ...
    'Left',left(:,order),'ConditionNumber',conditionNumber(order));
end

function [inputOperator,outputOperator,residual] = local_factor_l6_operator(b,n)
% Recover the per-pixel factorization B=R_out*L_L6.  For each L6 pixel,
% the corresponding S/C/I output rows must be collinear.
inputOperator=spalloc(n,3*n,ceil(nnz(b)/3));
outputOperator=spalloc(3*n,n,3*n);
for pixel=1:n
    rows=[pixel,n+pixel,2*n+pixel];
    block=b(rows,:);
    rowNorm=sqrt(full(sum(abs(block).^2,2)));
    [scale,reference]=max(rowNorm);
    if scale<=eps; continue; end
    direction=block(reference,:)/scale;
    coefficient=block*direction';
    inputOperator(pixel,:)=direction;
    outputOperator(rows,pixel)=coefficient;
end
residual=norm(b-outputOperator*inputOperator,'fro')/max(norm(b,'fro'),eps);
end

function local_plot_l6_factor_density(targetDensity,sourceDensity,totalDensity,data,outputDir)
maps=[targetDensity,sourceDensity,{totalDensity}];
labels={'target S','target C','target I','source S','source C','total'};
limit=max(cellfun(@(x)max(abs(x(:))),maps));
fig=figure('Visible','off','Color','w','Position',[80 80 1450 760]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for index=1:numel(maps)
    nexttile; imagesc(reshape(maps{index},data.MapSize)); axis image off; colorbar;
    if limit>0; clim([-limit limit]); end
    title(labels{index});
end
sgtitle('Physical L6-pixel attribution of Re(l^*Br)');
l6ns_save_figure(fig,outputDir,'physical_l6_pixel_attribution'); close(fig);
end

function local_plot_cross_actions(qA,qB,a,b,data,outputDir)
vectors={a*qA,b*qA,b*qB,a*qB};
labels={'Aq_A','Bq_A','Bq_B','Aq_B'};
populationLabels={'S','C','I','E'};
n=data.PopulationSize;
fig=figure('Visible','off','Color','w','Position',[60 60 1450 1250]);
tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
for row=1:4
    vector=vectors{row};
    maps={reshape(real(vector(1:n)),data.MapSize),reshape(real(vector(n+(1:n))),data.MapSize), ...
        reshape(real(vector(2*n+(1:n))),data.MapSize)};
    maps{4}=(1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit=max(cellfun(@(x)max(abs(x(:))),maps));
    for column=1:4
        nexttile; imagesc(maps{column}); axis image off; colorbar;
        if limit>0; clim([-limit limit]); end
        if row==1; title(populationLabels{column}); end
        if column==1; ylabel(labels{row},'Interpreter','none'); end
    end
end
sgtitle('How A and B act on each other''s leading positive modes');
l6ns_save_figure(fig,outputDir,'component_cross_actions'); close(fig);
end

function [atlas,summary,objects] = local_component_singular_analysis(operators,matrices,cfg,data,outputDir)
rows={}; summaryRows={}; objects=struct();
for operatorIndex=1:numel(operators)
    name=operators{operatorIndex}; matrix=matrices{operatorIndex};
    count=min(cfg.ComponentAtlasSize,size(matrix,1)-2);
    [uTop,sTop,vTop]=svds(matrix,count,'largest');
    sigmaTop=diag(sTop);
    nullAvailable=true;
    try
        [uNull,sNull,vNull]=svds(matrix,min(count,8),'smallest');
        sigmaNull=diag(sNull);
        [sigmaNull,order]=sort(sigmaNull,'ascend');
        uNull=uNull(:,order); vNull=vNull(:,order);
    catch exception
        warning('Near-null singular vectors unavailable for %s: %s',name,exception.message);
        uNull=[]; vNull=[]; sigmaNull=[]; nullAvailable=false;
    end
    frobenius=norm(matrix,'fro');
    stableRank=frobenius^2/max(sigmaTop(1)^2,eps);
    summaryRows(end+1,:)={name,sigmaTop(1),stableRank, ...
        sum(sigmaTop>sigmaTop(1)*1e-6),sum(sigmaTop>sigmaTop(1)*1e-8), ...
        min(sigmaTop),local_first_or_nan(sigmaNull),nullAvailable}; %#ok<AGROW>
    for modeIndex=1:numel(sigmaTop)
        inputFingerprint=local_fingerprint(vTop(:,modeIndex),data);
        outputFingerprint=local_fingerprint(uTop(:,modeIndex),data);
        rows(end+1,:)=local_singular_row(name,'leading','input',modeIndex,sigmaTop(modeIndex),inputFingerprint); %#ok<AGROW>
        rows(end+1,:)=local_singular_row(name,'leading','output',modeIndex,sigmaTop(modeIndex),outputFingerprint); %#ok<AGROW>
    end
    for modeIndex=1:numel(sigmaNull)
        inputFingerprint=local_fingerprint(vNull(:,modeIndex),data);
        outputFingerprint=local_fingerprint(uNull(:,modeIndex),data);
        rows(end+1,:)=local_singular_row(name,'near_null','input',modeIndex,sigmaNull(modeIndex),inputFingerprint); %#ok<AGROW>
        rows(end+1,:)=local_singular_row(name,'near_null','output',modeIndex,sigmaNull(modeIndex),outputFingerprint); %#ok<AGROW>
    end
    objects.(name)=struct('LeadingU',uTop,'LeadingV',vTop,'LeadingSigma',sigmaTop, ...
        'NullU',uNull,'NullV',vNull,'NullSigma',sigmaNull);
    local_plot_component_singular(uTop(:,1),vTop(:,1),uNull,vNull,data,name,outputDir);
end
atlas=cell2table(rows,'VariableNames',{'operator','category','side','mode','singularValue', ...
    'S_energyFraction','C_energyFraction','I_energyFraction','E_energyFraction', ...
    'participationRatio','dominantSpatialFrequency','dominantSpatialOrientationDeg'});
summary=cell2table(summaryRows,'VariableNames',{'operator','largestSingularValue','stableRank', ...
    'topKCountRelative1e6','topKCountRelative1e8','smallestOfComputedLeading', ...
    'smallestComputedSingularValue','nearNullAvailable'});
end

function row=local_singular_row(operator,category,side,index,sigma,fingerprint)
row={operator,category,side,index,sigma,fingerprint.SFraction,fingerprint.CFraction, ...
    fingerprint.IFraction,fingerprint.EFraction,fingerprint.ParticipationRatio, ...
    fingerprint.DominantFrequency,fingerprint.DominantOrientation};
end

function value=local_first_or_nan(values)
if isempty(values); value=NaN; else; value=values(1); end
end

function local_plot_component_singular(u,v,uNull,vNull,data,name,outputDir)
n=data.PopulationSize;
vectors={v,u}; labels={'leading input','leading output'};
populationLabels={'S','C','I','E'};
if ~isempty(vNull)
    vectors{end+1}=sqrt(sum(abs(vNull).^2,2)); labels{end+1}='near-null input leverage';
    vectors{end+1}=sqrt(sum(abs(uNull).^2,2)); labels{end+1}='near-null output leverage';
end
fig=figure('Visible','off','Color','w','Position',[60 60 1450 320*numel(vectors)]);
tiledlayout(numel(vectors),4,'TileSpacing','compact','Padding','compact');
for rowIndex=1:numel(vectors)
    vector=vectors{rowIndex};
    maps={reshape(real(vector(1:n)),data.MapSize),reshape(real(vector(n+(1:n))),data.MapSize), ...
        reshape(real(vector(2*n+(1:n))),data.MapSize)};
    maps{4}=(1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit=max(cellfun(@(x)max(abs(x(:))),maps));
    for column=1:4
        nexttile; imagesc(maps{column}); axis image off; colorbar;
        if limit>0 && rowIndex<=2; clim([-limit limit]); end
        if rowIndex==1; title(populationLabels{column}); end
        if column==1; ylabel(labels{rowIndex},'Interpreter','none'); end
    end
end
sgtitle(sprintf('%s singular input-output and near-null sectors',strrep(name,'_',' ')));
l6ns_save_figure(fig,outputDir,sprintf('singular_atlas_%s',name)); close(fig);
end

function fingerprint = local_fingerprint(vector,data)
n = data.PopulationSize;
s = vector(1:n); c = vector(n+(1:n)); i = vector(2*n+(1:n));
energy = [sum(abs(s).^2),sum(abs(c).^2),sum(abs(i).^2)];
total = sum(energy);
e = (1-data.ExcitatoryCWeight)*s+data.ExcitatoryCWeight*c;
eEnergy = sum(abs(e).^2);
pixelEnergy = abs(s).^2+abs(c).^2+abs(i).^2;
participation = sum(pixelEnergy)^2/max(sum(pixelEnergy.^2),eps);
[frequency,orientation] = local_dominant_fourier(reshape(e,data.MapSize));
fingerprint = struct('SFraction',energy(1)/total,'CFraction',energy(2)/total, ...
    'IFraction',energy(3)/total,'EFraction',eEnergy/(eEnergy+energy(3)), ...
    'ParticipationRatio',participation,'DominantFrequency',frequency, ...
    'DominantOrientation',orientation);
end

function [frequency,orientation] = local_dominant_fourier(map)
power = abs(fftshift(fft2(map))).^2;
[ny,nx] = size(power);
cy = floor(ny/2)+1; cx = floor(nx/2)+1;
power(cy,cx)=0;
[~,index]=max(power(:));
[iy,ix]=ind2sub(size(power),index);
kx=(ix-cx)/nx; ky=(iy-cy)/ny;
frequency=hypot(kx,ky);
orientation=mod(atan2d(ky,kx),180);
end

function tableOut = local_subspace_comparison(qA,qB,qJ)
labels={'A_to_B','A_to_J','B_to_J'};
pairs={qA,qB;qA,qJ;qB,qJ};
rows=cell(3,1);
for i=1:3
    cosine=svd(pairs{i,1}'*pairs{i,2});
    rows{i}={labels{i},median(cosine),min(cosine),max(cosine), ...
        sqrt(max(0,size(pairs{i,1},2)-norm(pairs{i,1}'*pairs{i,2},'fro')^2))};
end
tableOut=cell2table(vertcat(rows{:}),'VariableNames', ...
    {'comparison','medianCosine','minimumCosine','maximumCosine','grassmannDistance'});
end

function tableOut = local_redirection(operator,basis,label)
rows=cell(size(basis,2),1);
for i=1:size(basis,2)
    action=operator*basis(:,i);
    inside=basis*(basis'*action);
    outside=action-inside;
    rows{i}={label,i,norm(action),norm(inside),norm(outside), ...
        norm(outside)/max(norm(action),eps)};
end
tableOut=cell2table(vertcat(rows{:}),'VariableNames', ...
    {'direction','mode','actionNorm','insideNorm','outsideNorm','outsideFraction'});
end

function local_plot_mode_atlas(vectors,lambda,data,operatorName,category,outputDir)
n=data.PopulationSize;
fig=figure('Visible','off','Color','w','Position',[60 60 1400 1050]);
populationLabels={'S','C','I','E'};
tiledlayout(size(vectors,2),4,'TileSpacing','compact','Padding','compact');
for row=1:size(vectors,2)
    vector=vectors(:,row);
    maps={reshape(real(vector(1:n)),data.MapSize),reshape(real(vector(n+(1:n))),data.MapSize), ...
        reshape(real(vector(2*n+(1:n))),data.MapSize)};
    maps{4}=(1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit=max(cellfun(@(x)max(abs(x(:))),maps));
    for col=1:4
        nexttile; imagesc(maps{col}); axis image off; colorbar;
        if limit>0; clim([-limit limit]); end
        if row==1; title(populationLabels{col}); end
        if col==1; ylabel(sprintf('%d: %.3g%+.3gi',row,real(lambda(row)),imag(lambda(row)))); end
    end
end
sgtitle(sprintf('%s %s modes',strrep(operatorName,'_',' '),strrep(category,'_',' ')));
l6ns_save_figure(fig,outputDir,sprintf('atlas_%s_%s',operatorName,category)); close(fig);
end

function local_plot_attribution_scatter(t,outputDir)
fig=figure('Visible','off','Color','w','Position',[100 100 850 620]);
finiteCondition=t.conditionNumber(isfinite(t.conditionNumber) & t.conditionNumber>0);
if isempty(finiteCondition); finiteCondition=1; end
conditionCap=prctile(finiteCondition,99);
conditionForPlot=min(max(t.conditionNumber,1),conditionCap);
conditionForPlot(~isfinite(conditionForPlot))=conditionCap;
sizeValue=20+80*log10(conditionForPlot)./max(log10(conditionCap),1);
scatter(t.baseRadial,t.l6Radial,sizeValue,t.lambdaReal,'filled'); colorbar;
xline(0,'k--'); yline(0,'k--'); grid on;
xlabel('base radial contribution'); ylabel('L6 radial contribution');
title('Mode-by-mode stabilizing and destabilizing actions');
l6ns_save_figure(fig,outputDir,'full_mode_component_attribution'); close(fig);
end

function local_plot_complex_mode(vector,lambda,data,operatorName,category,outputDir)
n=data.PopulationSize;
population={vector(1:n),vector(n+(1:n)),vector(2*n+(1:n))};
population{4}=(1-data.ExcitatoryCWeight)*population{1}+data.ExcitatoryCWeight*population{2};
populationLabels={'S','C','I','E'};
fig=figure('Visible','off','Color','w','Position',[60 60 1450 1250]);
tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
rowLabels={'real','imaginary','magnitude','phase'};
for row=1:4
    for column=1:4
        map=reshape(population{column},data.MapSize);
        switch row
            case 1; shown=real(map);
            case 2; shown=imag(map);
            case 3; shown=abs(map);
            case 4; shown=angle(map);
        end
        nexttile; imagesc(shown); axis image off; colorbar;
        if row<=2
            limit=max(abs(shown(:))); if limit>0; clim([-limit limit]); end
        end
        if row==1; title(populationLabels{column}); end
        if column==1; ylabel(rowLabels{row}); end
    end
end
sgtitle(sprintf('%s %s complex mode, lambda %.4g%+.4gi', ...
    strrep(operatorName,'_',' '),strrep(category,'_',' '),real(lambda),imag(lambda)));
l6ns_save_figure(fig,outputDir,sprintf('complex_mode_%s_%s',operatorName,category)); close(fig);
end

function local_plot_redirection(t,outputDir)
fig=figure('Visible','off','Color','w','Position',[100 100 900 450]);
groups=unique(t.direction,'stable'); hold on;
for i=1:numel(groups)
    one=t(strcmp(t.direction,groups{i}),:);
    plot(one.mode,one.outsideFraction,'-o','LineWidth',1.4);
end
ylim([0 1]); grid on; xlabel('mode index'); ylabel('outside-subspace fraction');
legend(groups,'Interpreter','none','Location','best'); title('Component-induced mode rotation');
l6ns_save_figure(fig,outputDir,'component_mode_redirection'); close(fig);
end

function local_plot_adjoint_density(targetDensity,sourceDensity,data,outputDir)
n=data.PopulationSize;
fig=figure('Visible','off','Color','w','Position',[80 80 1450 720]);
populationLabels={'S','C','I','E'};
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for row=1:2
    if row==1
        density=targetDensity; label='target';
    else
        density=sourceDensity; label='source';
    end
    maps={reshape(density(1:n),data.MapSize),reshape(density(n+(1:n)),data.MapSize), ...
        reshape(density(2*n+(1:n)),data.MapSize)};
    maps{4}=(1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit=max(cellfun(@(x)max(abs(x(:))),maps));
    for col=1:4
        nexttile; imagesc(maps{col}); axis image off; colorbar;
        if limit>0; clim([-limit limit]); end
        title(sprintf('%s %s',label,populationLabels{col}));
    end
end
sgtitle('Exact adjoint localization of Re(l^* B r)');
l6ns_save_figure(fig,outputDir,'adjoint_l6_sensitivity_density'); close(fig);
end
