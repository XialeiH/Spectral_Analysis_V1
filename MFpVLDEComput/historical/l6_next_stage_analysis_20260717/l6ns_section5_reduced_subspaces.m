function result = l6ns_section5_reduced_subspaces(cfg,data,section1,section2)
% Section 5: useful shared subspaces and a validated reduced matrix pencil.

sectionDir=fullfile(cfg.OutputRoot,'section5_reduced_subspaces');
if ~exist(sectionDir,'dir'); mkdir(sectionDir); end

% Joint active directions from a normalized stacked operator.
gamma=norm(data.A,'fro')/max(norm(data.B,'fro'),eps);
stacked=[data.A;gamma*data.B];
jointCount=max(cfg.ReducedDimensions);
[~,jointS,jointV]=svds(stacked,jointCount,'largest');
jointSigma=diag(jointS);
activityRows=cell(jointCount,1);
for i=1:jointCount
    q=jointV(:,i);
    actionA=norm(data.A*q);
    actionB=gamma*norm(data.B*q);
    label="jointly_active";
    ratio=actionA/max(actionB,eps);
    if ratio>3; label="base_dominated";
    elseif ratio<1/3; label="L6_dominated";
    end
    activityRows{i}={'leading',i,jointSigma(i),actionA,actionB,ratio,label};
end
try
    [~,jointNullS,jointNullV]=svds(stacked,min(8,size(stacked,2)-2),'smallest');
    jointNullSigma=diag(jointNullS);
    [jointNullSigma,order]=sort(jointNullSigma,'ascend');
    jointNullV=jointNullV(:,order);
    for i=1:numel(jointNullSigma)
        q=jointNullV(:,i);
        actionA=norm(data.A*q); actionB=gamma*norm(data.B*q);
        ratio=actionA/max(actionB,eps);
        activityRows{end+1}={'near_null',i,jointNullSigma(i),actionA,actionB,ratio,'nearly_null_both'}; %#ok<AGROW>
    end
catch exception
    warning('L6NS:JointNullSVD','Joint near-null SVD unavailable: %s',exception.message);
end
activityTable=cell2table(vertcat(activityRows{:}),'VariableNames', ...
    {'sector','mode','stackedSingularValue','baseAction','normalizedL6Action','baseToL6Ratio','class'});
writetable(activityTable,fullfile(sectionDir,'joint_active_subspace.tsv'), ...
    'FileType','text','Delimiter','\t');
local_plot_joint_active(jointV(:,1:8),activityTable(1:8,:),data,sectionDir);

% Joint Krylov closure seeded by the persistent critical subspace and B edge.
qConsensus=section2.ConsensusObjects.k8.Basis;
bModes=l6ns_eigenpairs(data.B,8,'largestreal',cfg);
maxDimension=max(cfg.ReducedDimensions);
criticalDirection=section1.CriticalRight/norm(section1.CriticalRight);
remainingSeed=[qConsensus,bModes.Right,jointV(:,1:min(16,size(jointV,2)))];
remainingSeed=remainingSeed-criticalDirection*(criticalDirection'*remainingSeed);
[seedU,seedS,~]=svd(remainingSeed,'econ');
seedSingular=diag(seedS);
seedRank=sum(seedSingular>seedSingular(1)*1e-10);
q=[criticalDirection,seedU(:,1:min(seedRank,maxDimension-1))];
for iteration=1:8
    if size(q,2)>=maxDimension; break; end
    candidate=[data.A*q,data.B*q];
    candidate=candidate-q*(q'*candidate);
    [u,s,~]=svd(candidate,'econ');
    singular=diag(s);
    retained=sum(singular>singular(1)*1e-10);
    added=min([retained,maxDimension-size(q,2)]);
    if added==0; break; end
    q=[q,u(:,1:added)]; %#ok<AGROW>
end
closureBasis=q;

closureRows=cell(numel(cfg.ReducedDimensions),1);
for di=1:numel(cfg.ReducedDimensions)
    k=min(cfg.ReducedDimensions(di),size(closureBasis,2));
    qk=closureBasis(:,1:k);
    [leakA,leakB,fValue]=local_joint_leakage(qk,data.A,data.B);
    commutatorAction=data.A*(data.B*qk)-data.B*(data.A*qk);
    commutatorInside=qk' * commutatorAction;
    commutatorOutside=commutatorAction-qk*commutatorInside;
    denominator=max(norm(commutatorAction,'fro'),eps);
    closureRows{di}={k,leakA,leakB,fValue,norm(commutatorInside,'fro')/denominator, ...
        norm(commutatorOutside,'fro')/denominator};
end
closureTable=cell2table(vertcat(closureRows{:}),'VariableNames', ...
    {'dimension','leakageA','leakageB','jointObjective', ...
    'commutatorInsideFraction','commutatorOutsideFraction'});
writetable(closureTable,fullfile(sectionDir,'joint_krylov_closure.tsv'), ...
    'FileType','text','Delimiter','\t');
local_plot_closure(closureTable,sectionDir);

% Ordered critical spectral basis. The joint-Krylov basis above diagnoses
% closure, but its Galerkin compression can have spurious Ritz values for a
% strongly nonnormal operator. This nested basis is invariant for J(w_c),
% so every requested truncation contains the corresponding leading modes.
criticalMatrix=data.A+(1-section1.Summary.directWc)*data.B;
criticalModes=l6ns_eigenpairs(criticalMatrix,maxDimension,'largestreal',cfg);
[criticalBasis,~]=qr(criticalModes.Right,0);
criticalBasisRows=cell(numel(cfg.ReducedDimensions),1);
for di=1:numel(cfg.ReducedDimensions)
    k=min(cfg.ReducedDimensions(di),size(criticalBasis,2));
    qk=criticalBasis(:,1:k);
    reducedCritical=qk'*criticalMatrix*qk;
    invariantResidual=norm(criticalMatrix*qk-qk*reducedCritical,'fro')/ ...
        max(norm(criticalMatrix*qk,'fro'),eps);
    criticalBasisRows{di}={k,invariantResidual,max(real(eig(reducedCritical,'vector')))};
end
criticalBasisTable=cell2table(vertcat(criticalBasisRows{:}),'VariableNames', ...
    {'dimension','invariantResidualAtCritical','reducedMaxRealAtCritical'});
writetable(criticalBasisTable,fullfile(sectionDir,'critical_spectral_basis.tsv'), ...
    'FileType','text','Delimiter','\t');

% L6 coupling in the biorthogonal base eigenbasis.
modalCount=cfg.BaseModalCount;
aModes=l6ns_eigenpairs(data.A,modalCount,'largestreal',cfg);
crossGram=aModes.Left'*aModes.Right;
dualLeft=aModes.Left/crossGram';
biorthogonalityError=norm(dualLeft'*aModes.Right-eye(modalCount),'fro');
modalCoupling=dualLeft' * data.B * aModes.Right;
diagonalEnergy=norm(diag(diag(modalCoupling)),'fro')^2;
totalEnergy=norm(modalCoupling,'fro')^2;
couplingSummary=table(modalCount,biorthogonalityError,diagonalEnergy/totalEnergy, ...
    1-diagonalEnergy/totalEnergy,max(abs(diag(modalCoupling))), ...
    max(abs(modalCoupling-diag(diag(modalCoupling))),[],'all'), ...
    'VariableNames',{'modeCount','biorthogonalityError','diagonalEnergyFraction','offDiagonalEnergyFraction', ...
    'maximumDiagonalMagnitude','maximumOffDiagonalMagnitude'});
writetable(couplingSummary,fullfile(sectionDir,'base_basis_l6_coupling_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
local_plot_modal_coupling(modalCoupling,aModes.Lambda,sectionDir);

% Reduced pencils Q^*AQ+(1-w)Q^*BQ and validation against the full family.
validationW=unique([cfg.SubspaceWGrid(:);section1.Summary.directWc]);
reducedRows={};
for di=1:numel(cfg.ReducedDimensions)
    k=min(cfg.ReducedDimensions(di),size(criticalBasis,2));
    qk=criticalBasis(:,1:k);
    ar=qk'*data.A*qk;
    br=qk'*data.B*qk;
    reducedWc=local_dense_boundary(ar,br,cfg.CriticalBoundary);
    for wi=1:numel(validationW)
        w=validationW(wi);
        jr=ar+(1-w)*br;
        [vr,dr]=eig(jr,'vector');
        [reducedMax,index]=max(real(dr));
        lifted=qk*vr(:,index); lifted=lifted/norm(lifted);
        fullModes=l6ns_eigenpairs(data.A+(1-w)*data.B,8,'largestreal',cfg);
        [fullMax,fullIndex]=max(real(fullModes.Lambda));
        overlap=abs(fullModes.Right(:,fullIndex)'*lifted)/ ...
            max(norm(fullModes.Right(:,fullIndex))*norm(lifted),eps);
        liftedAll=qk*vr;
        modeOverlaps=abs(fullModes.Right(:,fullIndex)'*liftedAll) ./ ...
            max(norm(fullModes.Right(:,fullIndex))*vecnorm(liftedAll),eps);
        [targetOverlap,targetIndex]=max(modeOverlaps);
        targetLambda=dr(targetIndex);
        reducedRows(end+1,:)={k,w,fullMax,reducedMax,abs(fullMax-reducedMax), ...
            overlap,real(targetLambda),imag(targetLambda), ...
            abs(targetLambda-fullModes.Lambda(fullIndex)),targetOverlap, ...
            reducedWc,reducedWc-section1.Summary.directWc}; %#ok<AGROW>
    end
end
reducedTable=cell2table(reducedRows,'VariableNames', ...
    {'dimension','w','fullMaxReal','reducedMaxReal','absoluteLeadingError', ...
    'liftedLeadingModeOverlap','targetModeReal','targetModeImag', ...
    'targetModeAbsoluteError','targetModeOverlap','reducedWCritical','wCriticalError'});
writetable(reducedTable,fullfile(sectionDir,'reduced_pencil_validation.tsv'), ...
    'FileType','text','Delimiter','\t');
local_plot_reduced_validation(reducedTable,section1.Summary.directWc,sectionDir);

% Validate continuous-time transient gain and lifted input-output patterns.
transientValidation=local_reduced_transient_validation(cfg,data,section1,section2,criticalBasis);
writetable(transientValidation,fullfile(sectionDir,'reduced_transient_validation.tsv'), ...
    'FileType','text','Delimiter','\t');

result=struct('JointActive',activityTable,'JointActiveBasis',jointV, ...
    'Closure',closureTable,'ClosureBasis',closureBasis,'ModalCoupling',modalCoupling, ...
    'CriticalBasisSummary',criticalBasisTable,'CriticalBasis',criticalBasis, ...
    'CouplingSummary',couplingSummary,'ReducedValidation',reducedTable, ...
    'TransientValidation',transientValidation);
save(fullfile(sectionDir,'section5_result.mat'),'result','-v7.3');
end

function [leakA,leakB,fValue]=local_joint_leakage(q,a,b)
actionA=a*q; actionB=b*q;
outsideA=actionA-q*(q'*actionA);
outsideB=actionB-q*(q'*actionB);
leakA=norm(outsideA,'fro')/max(norm(actionA,'fro'),eps);
leakB=norm(outsideB,'fro')/max(norm(actionB,'fro'),eps);
fValue=leakA^2+leakB^2;
end

function wc=local_dense_boundary(a,b,boundary)
highW=0; lowW=-0.2;
f=@(w)max(real(eig(a+(1-w)*b)));
while f(lowW)<boundary && lowW>-20; lowW=1.75*lowW; end
if f(highW)>=boundary || f(lowW)<boundary
    wc=NaN;
    return
end
for iteration=1:45
    middle=0.5*(lowW+highW);
    if f(middle)>=boundary; lowW=middle; else; highW=middle; end
end
wc=0.5*(lowW+highW);
end

function tableOut=local_reduced_transient_validation(cfg,data,section1,section2,closureBasis)
wValues=unique([0;1;section1.Summary.directWc]);
rows={};
for dimension=cfg.ReducedDimensions
    k=min(dimension,size(closureBasis,2));
    q=closureBasis(:,1:k);
    for wi=1:numel(wValues)
        w=wValues(wi);
        field=local_w_field(w);
        fullTransient=section2.TransientObjects.(field);
        weights=fullTransient.WeightVector;
        qScaled=orth(weights.*q);
        generator=data.A+(1-w)*data.B-speye(data.Dimension);
        scaledAction=weights.*(generator*(qScaled./weights));
        reducedGenerator=qScaled'*scaledAction;
        gain=nan(numel(cfg.TransientTimeGrid),1);
        vectors=cell(numel(cfg.TransientTimeGrid),2);
        for ti=1:numel(cfg.TransientTimeGrid)
            propagator=expm(cfg.TransientTimeGrid(ti)*reducedGenerator);
            [u,s,v]=svd(propagator,'econ');
            gain(ti)=s(1,1);
            vectors{ti,1}=v(:,1); vectors{ti,2}=u(:,1);
        end
        [peakGain,peakIndex]=max(gain);
        liftedInput=(qScaled*vectors{peakIndex,1})./weights;
        liftedOutput=(qScaled*vectors{peakIndex,2})./weights;
        rows(end+1,:)={k,w,fullTransient.PeakGain,peakGain, ...
            abs(peakGain-fullTransient.PeakGain)/max(fullTransient.PeakGain,eps), ...
            fullTransient.PeakTime,cfg.TransientTimeGrid(peakIndex), ...
            abs((weights.*liftedInput)'*(weights.*fullTransient.OptimalInput))/ ...
            max(norm(weights.*liftedInput)*norm(weights.*fullTransient.OptimalInput),eps), ...
            abs((weights.*liftedOutput)'*(weights.*fullTransient.PeakOutput))/ ...
            max(norm(weights.*liftedOutput)*norm(weights.*fullTransient.PeakOutput),eps)}; %#ok<AGROW>
    end
end
tableOut=cell2table(rows,'VariableNames',{'dimension','w','fullPeakGain','reducedPeakGain', ...
    'relativePeakGainError','fullPeakTime','reducedPeakTime', ...
    'liftedInputOverlap','liftedOutputOverlap'});
end

function field=local_w_field(w)
field=sprintf('w_%+.6f',w);
field=strrep(field,'+','p'); field=strrep(field,'-','m'); field=strrep(field,'.','p');
end

function local_plot_joint_active(q,t,data,outputDir)
n=data.PopulationSize;
fig=figure('Visible','off','Color','w','Position',[60 60 1400 1050]);
populationLabels={'S','C','I','E'};
tiledlayout(size(q,2),4,'TileSpacing','compact','Padding','compact');
for row=1:size(q,2)
    maps={reshape(real(q(1:n,row)),data.MapSize),reshape(real(q(n+(1:n),row)),data.MapSize), ...
        reshape(real(q(2*n+(1:n),row)),data.MapSize)};
    maps{4}=(1-data.ExcitatoryCWeight)*maps{1}+data.ExcitatoryCWeight*maps{2};
    limit=max(cellfun(@(x)max(abs(x(:))),maps));
    for col=1:4
        nexttile; imagesc(maps{col}); axis image off; colorbar;
        if limit>0; clim([-limit limit]); end
        if row==1; title(populationLabels{col}); end
        if col==1; ylabel(sprintf('%d %s',row,t.class(row)),'Interpreter','none'); end
    end
end
sgtitle('Leading jointly active directions of [A; \gamma B]');
l6ns_save_figure(fig,outputDir,'joint_active_direction_maps'); close(fig);
end

function local_plot_closure(t,outputDir)
fig=figure('Visible','off','Color','w','Position',[100 100 1050 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(t.dimension,t.leakageA,'-o','LineWidth',1.4); hold on;
plot(t.dimension,t.leakageB,'-o','LineWidth',1.4); plot(t.dimension,t.jointObjective,'-o','LineWidth',1.4);
grid on; xlabel('dimension'); ylabel('leakage/objective'); legend({'A','B','F_k'},'Location','best');
title('Joint Krylov closure');
nexttile; plot(t.dimension,t.commutatorInsideFraction,'-o','LineWidth',1.4); hold on;
plot(t.dimension,t.commutatorOutsideFraction,'-o','LineWidth',1.4); grid on;
xlabel('dimension'); ylabel('commutator action fraction'); legend({'inside','outside'},'Location','best');
title('Relevant-subspace noncommutation');
l6ns_save_figure(fig,outputDir,'joint_closure_and_commutator'); close(fig);
end

function local_plot_modal_coupling(coupling,lambda,outputDir)
fig=figure('Visible','off','Color','w','Position',[100 100 1000 760]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; imagesc(log10(abs(coupling)+1e-14)); axis image; colorbar;
xlabel('base right mode'); ylabel('base left mode'); title('log_{10}|L_A^*BR_A|');
nexttile; plot(real(lambda),'-o'); grid on; xlabel('base mode index'); ylabel('Re \lambda(A)');
title('Ordering of the base modes');
l6ns_save_figure(fig,outputDir,'l6_coupling_in_base_basis'); close(fig);
end

function local_plot_reduced_validation(t,fullWc,outputDir)
fig=figure('Visible','off','Color','w','Position',[100 100 1100 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
dimensions=unique(t.dimension);
for i=1:numel(dimensions)
    one=t(t.dimension==dimensions(i),:);
    semilogy(one.w,one.absoluteLeadingError,'-o','LineWidth',1.3);
end
grid on; xlabel('w'); ylabel('|maxRe full - reduced|'); legend(string(dimensions),'Location','best');
title('Leading-eigenvalue accuracy');
nexttile;
    dimensions=unique(t.dimension);
    thresholdError=nan(size(dimensions));
    for i=1:numel(dimensions)
        first=find(t.dimension==dimensions(i),1,'first');
        thresholdError(i)=t.wCriticalError(first);
    end
    plot(dimensions,abs(thresholdError),'-o','LineWidth',1.5); hold on;
yline(0,'k--'); grid on; xlabel('dimension'); ylabel('|w_c^{red}-w_c|');
title(sprintf('Threshold accuracy, full w_c=%.6f',fullWc));
l6ns_save_figure(fig,outputDir,'reduced_pencil_validation'); close(fig);
end
