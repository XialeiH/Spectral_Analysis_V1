function compute_supplement_S1_S2(whichFigure)
% Recompute supplementary comparisons from the preserved h96 and genuine tables.
if nargin<1, whichFigure='both'; end
root='/Users/xialeihuang/Desktop/Neuroscience_Project';
model=fullfile(root,'Spectral_Analysis','matlab-inserting_into_CG_model');
repo=fullfile(model,'Complete_Code_for_Paper3','NYU-Vision-2Drive-main');
analysis=fullfile(model,'l6_next_stage_analysis_20260717');
out=fullfile(model,'report_figure_artifacts','Supplement_S1_S2');
if ~isfolder(out), mkdir(out); end
addpath(fullfile(repo,'Utils'));
addpath(analysis,'-begin');
runtime=fullfile(analysis,'real_tuning_l6_i_20260722','mechanism_analysis_20260722','runtime_h96');
addpath(runtime,'-begin');
assert(startsWith(which('predict_pref6D_S'),runtime));
loaded=load(fullfile(analysis,'results_global_bifurcation_20260722','global_bifurcation_setup.mat'),'setup');
context=loaded.setup.Context;
assert(numel(context.FixedPoint)==4800 && all(context.OrientationUse==0,'all') && all(context.ContrastUse==100,'all'));
dataRoot=fullfile(repo,'Data','Paper2_NetworkTuning','Fig1V4');
if any(strcmp(whichFigure,{'S1','both'}))
    tauTable=readtable(fullfile(analysis,'results_tau_match','tau_fit_lif_peak','tau_lif_fit_summary.tsv'), ...
        'FileType','text','Delimiter','\t');
    tauMs=tauTable.tauMs(1); lifPeakMs=tauTable.lifPeakTimeMs(1);
    initialState=context.InitialState(:);
    tables=load_tables(dataRoot,3);
    phiLibrary=@(x) library3(x,context,tables);
    % Check our explicit map against one step of the original library driver.
    original=LDEIteration_135FuncMain_CombDom_RealLGNL6(context.PixLGNCtgr, ...
        context.L6Kernel,split_state(initialState),1,context.L6Parameters,1, ...
        context.C_SS,context.C_CS,context.C_IS,context.C_SC,context.C_CC,context.C_IC, ...
        context.C_SI,context.C_CI,context.C_II,context.L4SEp,context.L4SIp, ...
        context.L4CEp,context.L4CIp,context.L4IEp,context.L4IIp, ...
        {tables.x},{tables.y},{tables.raw},(1:4)',4,10,10,context.Isaturation,'xn', ...
        context.EKpUse,context.IKpUse);
    reference=[original{2}.S;original{2}.C;original{2}.I];
    originalStepError=norm(reference-phiLibrary(initialState))/norm(reference);
    assert(originalStepError<1e-12);
    hValues=[0.02 0.01];
    iterations=cell(2,1);
    for k=1:2
        h=hValues(k); t=(0:ceil(180/(h*tauMs)))*h*tauMs;
        x=initialState; means=zeros(numel(t),4); pixels=zeros(numel(t),2);
        for j=1:numel(t)
            [means(j,:),pixels(j,:)]=observables(x);
            if j<numel(t), x=x+h*(phiLibrary(x)-x); end
        end
        iterations{k}=struct('h',h,'timeMs',t(:),'means',means,'pixels',pixels,'lastState',x);
        fprintf('S1 library h=%g finished: %d steps\n',h,numel(t)-1);
    end
    phiDnn=@(x) l6ns_phi(x,0,context);
    odeTimesMs=(0:0.1:180)';
    opts=odeset('RelTol',2e-7,'AbsTol',1e-9,'MaxStep',0.2);
    [~,states]=ode45(@(~,x) (phiDnn(x)-x)/tauMs,odeTimesMs,initialState,opts);
    odeMeans=zeros(numel(odeTimesMs),4); odePixels=zeros(numel(odeTimesMs),2);
    for j=1:numel(odeTimesMs), [odeMeans(j,:),odePixels(j,:)]=observables(states(j,:)'); end
    library=iterations{2};
    [libraryPeakMs,libraryPeakRate]=peak_time(library.timeMs,library.means(:,4));
    [odePeakMs,odePeakRate]=peak_time(odeTimesMs,odeMeans(:,4));
    coarse=iterations{1};
    coarseOnFine=interp1(coarse.timeMs,coarse.means,library.timeMs,'linear');
    stepRefinementMaxRateError=max(abs(coarseOnFine-library.means),[],'all');
    summary=table(tauMs,lifPeakMs,library.h,library.h*tauMs,libraryPeakMs,odePeakMs, ...
        libraryPeakRate,odePeakRate,stepRefinementMaxRateError,originalStepError, ...
        'VariableNames',{'tauMs','lifPeakMs','h','iterationStepMs','libraryPeakMs','odePeakMs', ...
        'libraryPeakMeanE','odePeakMeanE','stepRefinementMaxRateError','originalStepError'});
    disp(summary); writetable(summary,fullfile(out,'S1_summary.tsv'),'FileType','text','Delimiter','\t');
    save(fullfile(out,'S1_data.mat'),'library','coarse','odeTimesMs','odeMeans','odePixels', ...
        'initialState','summary','tauMs','lifPeakMs','-v7.3');
end
if any(strcmp(whichFigure,{'S1','S1_original','both'}))
    saved=load(fullfile(out,'S1_data.mat'),'initialState','tauMs','summary');
    if ~exist('tables','var'), tables=load_tables(dataRoot,3); end
    h=0.33; t=(0:ceil(180/(h*saved.tauMs)))*h*saved.tauMs;
    x=saved.initialState; means=zeros(numel(t),4); pixels=zeros(numel(t),2);
    for j=1:numel(t)
        [means(j,:),pixels(j,:)]=observables(x);
        if j<numel(t), x=x+h*(library3(x,context,tables)-x); end
    end
    originalLibrary=struct('h',h,'timeMs',t(:),'means',means,'pixels',pixels,'lastState',x);
    [peakMs,peakRate]=peak_time(t(:),means(:,4));
    summary=saved.summary;
    summary.originalLibraryPeakMs=peakMs; summary.originalLibraryPeakMeanE=peakRate;
    summary.originalLibraryStepMs=h*saved.tauMs;
    save(fullfile(out,'S1_data.mat'),'originalLibrary','summary','-append');
    writetable(summary,fullfile(out,'S1_summary.tsv'),'FileType','text','Delimiter','\t');
    disp(summary);
end
if any(strcmp(whichFigure,{'S2_moved','S2','both'}))
    phiMoved=@(x) l6ns_phi(x,1,context,[1 1],1,'true');
    fixed=context.FixedPoint(:); h=0.1;
    for iteration=1:20000
        target=phiMoved(fixed);
        residual=norm(target-fixed)/max(norm(fixed),eps);
        if residual<1e-10, break; end
        fixed=fixed+h*(target-fixed);
    end
    assert(residual<1e-8,'DNN zero-L6 fixed point did not converge.');
    zeroSetup=loaded.setup;
    zeroSetup.Context.FixedL6LibInd=zeros(size(context.FixedL6LibInd));
    addpath(fullfile(model,'pixel_size_scaling_20260813','runtime_override'),'-end');
    J=l6ns_state_jacobian(zeroSetup,'L6',fixed,1);
    rng(2410); direction=randn(4800,1); direction=direction/norm(direction);
    derivativeErrors=zeros(3,1);
    for k=1:3
        step=1e-3/2^(k-1);
        fd=(phiMoved(fixed+step*direction)-phiMoved(fixed-step*direction))/(2*step);
        derivativeErrors(k)=norm(J*direction-fd)/max(norm(fd),eps);
    end
    assert(max(derivativeErrors)<1e-4);
    moved=mode_summary(J,fixed);
    moved.fixedPointResidual=residual; moved.derivativeErrors=derivativeErrors;
    moved.iterationCount=iteration; moved.iterationH=h; moved.fixedL6Input=0;
    moved.fixedPointShift=norm(fixed-context.FixedPoint(:));
    save(fullfile(out,'S2_dnn_moved.mat'),'moved','J','-v7.3');
    fprintf('S2 DNN zero-L6 moved: maxReal=%g singular=%g residual=%g derivative error=%g steps=%d\n', ...
        max(real(moved.eigenvalues)),moved.singularValue,residual,max(derivativeErrors),iteration);
end
if any(strcmp(whichFigure,{'S2','both'}))
    % FPP removes only the L6 derivative, not its fixed-point input value.
    J=loaded.setup.Pathway.JBaseline-loaded.setup.Pathway.J6;
    fixed=context.FixedPoint(:);
    phiFrozen=@(x) l6ns_phi(x,1,context);
    residual=norm(phiFrozen(fixed)-fixed)/norm(fixed);
    assert(residual<1e-9);
    rng(2409); direction=randn(4800,1); direction=direction/norm(direction);
    step=1e-3; fd=(phiFrozen(fixed+step*direction)-phiFrozen(fixed-step*direction))/(2*step);
    derivativeError=norm(J*direction-fd)/norm(fd); assert(derivativeError<1e-4);
    dnn=mode_summary(J,fixed);
    dnn.fixedPointResidual=residual; dnn.derivativeError=derivativeError;
    save(fullfile(out,'S2_dnn.mat'),'dnn','J','-v7.3');
    fprintf('S2 DNN frozen: maxReal=%g residual=%g derivative error=%g\n',max(real(dnn.eigenvalues)),residual,derivativeError);
    clear J loaded
    p=load(fullfile(dataRoot,'AllMFPixPara_Paper2TuneFig1V4D2.mat'), ...
        'N_HC','n_S_HC','NnSPixel','OD_SMap','C_SS_Pixel_Us','C_CS_Pixel_Us', ...
        'C_IS_Pixel_Us','C_SC_Pixel_Us','C_CC_Pixel_Us','C_IC_Pixel_Us', ...
        'C_SI_Pixel_Us','C_CI_Pixel_Us','C_II_Pixel_Us');
    pairs={'SS','CS','IS','SC','CC','IC','SI','CI','II'};
    for k=1:numel(pairs)
        cc.(['C_' pairs{k}])=sparse(AveSpatKer(p.(['C_' pairs{k} '_Pixel_Us']),p.N_HC,4,10,10));
    end
    l6Filter=SpatialGaussianFilt_my(p.OD_SMap,p.N_HC,p.n_S_HC,p.n_S_HC*0.34,1.25,false);
    lgnFilter=SpatialGaussianFilt_my(p.OD_SMap,p.N_HC,p.n_S_HC,p.n_S_HC*0.2,1,false);
    l6Weights=LGNIndSpat(l6Filter,1:4,p.NnSPixel,p.N_HC,4,10,10);
    lgnWeights=LGNIndSpat(lgnFilter,1:4,p.NnSPixel,p.N_HC,4,10,10);
    weights=reshape(lgnWeights,1600,4,1).*reshape(l6Weights,1600,1,4);
    tables=load_tables(dataRoot,2);
    cc.AS=cc.C_SS+cc.C_CS+cc.C_IS;
    cc.AC=cc.C_SC+cc.C_CC+cc.C_IC;
    cc.AI=cc.C_SI+cc.C_CI+cc.C_II;
    phi=@(x) library2(x,cc,tables,weights);
    x0=[2.5*ones(1600,1);8*ones(1600,1);16*ones(1600,1)];
    original=LDEIteration_16FuncMain_CombDom(weights,split_state(x0),1,1, ...
        cc.C_SS,cc.C_CS,cc.C_IS,cc.C_SC,cc.C_CC,cc.C_IC,cc.C_SI,cc.C_CI,cc.C_II, ...
        {tables.x},{tables.y},{tables.raw},4,10,10,false,'xn');
    reference=[original{2}.S;original{2}.C;original{2}.I];
    originalStepError=norm(reference-phi(x0))/norm(reference); assert(originalStepError<1e-12);
    x=x0; h=0.1;
    for iteration=1:20000
        target=phi(x); residual=norm(target-x)/max(norm(x),eps);
        if residual<1e-10, break; end
        x=x+h*(target-x);
    end
    assert(residual<1e-8,'Paper2 fixed point did not converge.');
    fprintf('S2 Paper2 equilibrium: %d steps, residual=%g\n',iteration,residual);
    % Central finite differences of both library inputs, followed by the
    % exact linear connectivity chain. No surrogate derivatives are used.
    X=cc.AS*x(1:1600)+cc.AC*x(1601:3200); Y=cc.AI*x(3201:4800);
    epsInput=1e-3;
    gX=(response2(X+epsInput,Y,tables,weights)-response2(X-epsInput,Y,tables,weights))/(2*epsInput);
    gY=(response2(X,Y+epsInput,tables,weights)-response2(X,Y-epsInput,tables,weights))/(2*epsInput);
    blocks=cell(3,3);
    for k=1:3
        rows=(k-1)*1600+(1:1600);
        dx=spdiags(gX(rows),0,1600,1600); dy=spdiags(gY(rows),0,1600,1600);
        blocks(k,:)={dx*cc.AS,dx*cc.AC,dy*cc.AI};
    end
    J=cell2mat(blocks);
    derivativeErrors=zeros(3,1);
    for k=1:3
        step=1e-3/2^(k-1);
        fd=(phi(x+step*direction)-phi(x-step*direction))/(2*step);
        derivativeErrors(k)=norm(J*direction-fd)/norm(fd);
    end
    assert(max(derivativeErrors)<1e-4);
    paper2=mode_summary(J,x);
    paper2.fixedPointResidual=residual; paper2.derivativeErrors=derivativeErrors;
    paper2.originalStepError=originalStepError; paper2.iterationCount=iteration;
    paper2.iterationH=h; paper2.finiteDifferenceInputStep=epsInput;
    save(fullfile(out,'S2_paper2.mat'),'paper2','J','weights','cc','-v7.3');
    fprintf('S2 Paper2 maxReal=%g singular=%g numerical error=%g\n', ...
        max(real(paper2.eigenvalues)),paper2.singularValue,max(derivativeErrors));
end
end

function tables=load_tables(dataRoot,dimension)
domains={'Small','Large','LARGER'};
for k=1:3
    if dimension==3
        file=fullfile(dataRoot,'Paper3PlotingData', ...
            sprintf('Func200V4D2-L6RealCtrlL4-Ang0.0-SF2.5TF10-%s.mat',domains{k}));
    else
        file=fullfile(dataRoot,'Paper2PlotingData',sprintf('Func16V4D2Ang0.0%s.mat',domains{k}));
    end
    a=load(file,'L4EmeshX','L4ImeshY','LDEFrfunc');
    tables(k).x=a.L4EmeshX; tables(k).y=a.L4ImeshY; tables(k).raw=a.LDEFrfunc;
    tables(k).file=file;
    if dimension==3
        for population='SCI'
            cellData=a.LDEFrfunc.(population); sample=cellData{1,1};
            packed=zeros(size(cellData,1),size(sample,1),size(sample,2),size(cellData,2));
            for i=1:size(cellData,1)
                for j=1:size(cellData,2), packed(i,:,:,j)=cellData{i,j}; end
            end
            tables(k).packed.(population)=packed;
            for i=1:size(cellData,1)
                tables(k).interpolant.(population){i}=griddedInterpolant( ...
                    {unique(a.L4ImeshY(:)),unique(a.L4EmeshX(:)),(1:size(cellData,2))'}, ...
                    squeeze(packed(i,:,:,:)),'linear','none');
            end
        end
    end
end
end

function y=library3(x,c,tables)
s=x(1:1600); v=x(1601:3200); i=x(3201:4800);
if c.Isaturation
    e=.6923*s+.3077*v; scale=L6Convert(e,c.EKpUse)./e;
    scale(~isfinite(scale))=1; s=s.*scale; v=v.*scale; i=InhMulp(i,c.IKpUse);
end
X=[(c.C_SS*s+c.C_SC*v)/c.L4SEp,(c.C_CS*s+c.C_CC*v)/c.L4CEp,(c.C_IS*s+c.C_IC*v)/c.L4IEp];
Y=[c.C_SI*i/c.L4SIp,c.C_CI*i/c.L4CIp,c.C_II*i/c.L4IIp];
e=reshape(.6923*s+.3077*v,40,40);
filtered=conv2(padarray(e,[1 1],'circular'),c.L6Kernel,'same');
z=L6Convert(filtered(2:end-1,2:end-1),c.L6Parameters)/3;
n6=size(tables(1).raw.S,2); z=min(max(z(:),1),n6);
for d=1:numel(tables)
    values=zeros(1600,3); names='SCI';
    for k=1:3
        for category=1:size(c.PixLGNCtgr,2)
            values(:,k)=values(:,k)+c.PixLGNCtgr(:,category).* ...
                tables(d).interpolant.(names(k)){category}(Y(:,k),X(:,k),z);
        end
    end
    if all(isfinite(values),'all'), y=values(:); return; end
end
error('All three genuine Paper3 library domains were exhausted.');
end

function y=library2(x,c,tables,weights)
X=c.AS*x(1:1600)+c.AC*x(1601:3200); Y=c.AI*x(3201:4800);
y=response2(X,Y,tables,weights);
end

function y=response2(X,Y,tables,weights)
for d=1:numel(tables)
    values=zeros(1600,3); names='SCI';
    for k=1:3
        values(:,k)=LDEIterFunc_Grating_16Func(tables(d).x,tables(d).y, ...
            tables(d).raw.(names(k)),X,Y,weights);
    end
    if all(isfinite(values),'all'), y=values(:); return; end
end
error('All three genuine Paper2 library domains were exhausted.');
end

function result=mode_summary(J,fixed)
fprintf('Computing complete %d-state eigenspectrum...\n',size(J,1));
[vectors,values]=eig(full(J),'vector');
[~,order]=sortrows([real(values),imag(values)],[-1 -2]);
result.eigenvalues=values(order); v=vectors(:,order(1)); clear vectors
[~,pivot]=max(abs(v)); v=v*exp(-1i*angle(v(pivot))); v=v/norm(v);
opts=struct('tol',1e-10,'maxit',3000,'disp',0);
[u,s,r,flag]=svds(J,1,'largest',opts); assert(flag==0);
% Figure 1F's last column is the left/output singular mode U.
[~,pivot]=max(abs(u)); phase=exp(-1i*angle(u(pivot))); u=u*phase; r=r*phase;
result.fixedPoint=fixed; result.eigenvector=v;
result.leftSingularVector=u; result.rightSingularVector=r; result.singularValue=s;
result.fixedPointE=reshape(.6923*fixed(1:1600)+.3077*fixed(1601:3200),40,40);
result.eigenvectorE=reshape(real(.6923*v(1:1600)+.3077*v(1601:3200)),40,40);
result.singularModeE=reshape(real(.6923*u(1:1600)+.3077*u(1601:3200)),40,40);
result.eigenResidual=norm(J*v-result.eigenvalues(1)*v)/max(norm(J*v),eps);
result.singularResidual=norm(J*r-u*s)/max(s,eps);
assert(result.eigenResidual<1e-7 && result.singularResidual<1e-7);
end

function s=split_state(x)
s=struct('S',x(1:1600),'C',x(1601:3200),'I',x(3201:4800));
end

function [means,pixels]=observables(x)
s=split_state(x); e=.6923*s.S+.3077*s.C;
means=[mean(s.S),mean(s.C),mean(s.I),mean(e)];
pixels=e(sub2ind([40 40],[5 1],[10 10]))';
end

function [time,value]=peak_time(t,y)
[value,k]=max(y); time=t(k);
if k>1 && k<numel(y)
    q=polyfit(t(k-1:k+1)-t(k),y(k-1:k+1),2);
    offset=-q(2)/(2*q(1));
    if q(1)<0 && abs(offset)<=max(diff(t(k-1:k+1)))
        time=t(k)+offset; value=polyval(q,offset);
    end
end
end
