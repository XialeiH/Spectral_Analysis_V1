function run_followup_ode_hidden_mode()
% Integrate the nonlinear continuous-time ODE around the hidden L6 instability.

paths=followup_initialize();
loaded=load(paths.SetupFile,'setup');
setup=loaded.setup; context=setup.Context;
beta=0.17; gain6=1+beta;
fixed=followup_branch_state(paths,"L6",beta,"low",context);
previous=followup_branch_state(paths,"L6",beta-0.001,"low",context);
phi=@(x)mechanism_phi_variant(x,context,gain6,1,'extended');
fixed=local_refine(phi,fixed,context.RelaxationP);
J=real_tuning_true_jacobian(fixed,context,gain6,1);
[lambda,right,left,eigenResidual]=followup_leading_pair(J);

rng(2906,'twister');
directions=zeros(numel(fixed),3);
directions(:,1)=real(right)/norm(real(right));
directions(:,2)=(previous-fixed)/norm(previous-fixed);
directions(:,3)=randn(size(fixed));
directions(:,3)=directions(:,3)/norm(directions(:,3));
names=["leading full-space mode";"symmetric continuation direction"; ...
    "random symmetry-breaking direction"];
epsilon=1e-7*norm(fixed);
dtIntrinsic=0.25;
finalIntrinsic=550;
steps=round(finalIntrinsic/dtIntrinsic);
saveEvery=4;
saveSteps=(0:saveEvery:steps)';
intrinsicTime=saveSteps*dtIntrinsic;
timeMs=intrinsicTime*setup.Config.TauMs;
modalAmplitude=nan(numel(saveSteps),3);
relativeDistance=nan(numel(saveSteps),3);
initialProjection=nan(3,1);
for directionIndex=1:3
    state=fixed+epsilon*directions(:,directionIndex);
    initialProjection(directionIndex)=abs(left'*(state-fixed))/epsilon;
    outputIndex=1;
    modalAmplitude(outputIndex,directionIndex)= ...
        abs(left'*(state-fixed))/epsilon;
    relativeDistance(outputIndex,directionIndex)=norm(state-fixed)/epsilon;
    for step=1:steps
        k1=phi(state)-state;
        predictor=state+dtIntrinsic*k1;
        k2=phi(predictor)-predictor;
        state=state+(dtIntrinsic/2)*(k1+k2);
        if mod(step,saveEvery)==0
            outputIndex=outputIndex+1;
            modalAmplitude(outputIndex,directionIndex)= ...
                abs(left'*(state-fixed))/epsilon;
            relativeDistance(outputIndex,directionIndex)=norm(state-fixed)/epsilon;
        end
    end
end
linearPrediction=exp((real(lambda)-1)*intrinsicTime);

cells=cell(3,1);
for index=1:3
    cells{index}=table(repmat(index,numel(intrinsicTime),1), ...
        repmat(names(index),numel(intrinsicTime),1),intrinsicTime,timeMs, ...
        modalAmplitude(:,index),relativeDistance(:,index),linearPrediction, ...
        'VariableNames',{'directionId','direction','intrinsicTimeTauUnits', ...
        'timeMs','normalizedLeadingModeAmplitude','normalizedStateDistance', ...
        'linearPrediction'});
end
trajectoryTable=vertcat(cells{:});
summaryTable=table((1:3)',names,initialProjection,modalAmplitude(end,:)', ...
    relativeDistance(end,:)',repmat(real(lambda),3,1), ...
    repmat((real(lambda)-1)/setup.Config.TauMs,3,1), ...
    repmat(eigenResidual,3,1), ...
    'VariableNames',{'directionId','direction','initialLeadingModeProjection', ...
    'finalNormalizedLeadingModeAmplitude','finalNormalizedStateDistance', ...
    'leadingMapEigenvalue','predictedGrowthPerMs','eigenpairResidual'});
outDir=fullfile(paths.FollowupOutput,'ode_hidden_mode');
if ~exist(outDir,'dir'); mkdir(outDir); end
writetable(trajectoryTable,fullfile(outDir,'ode_hidden_mode_trajectory.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(summaryTable,fullfile(outDir,'ode_hidden_mode_summary.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outDir,'ode_hidden_mode.mat'),'trajectoryTable','summaryTable', ...
    'fixed','lambda','right','left','directions','-v7.3');
fprintf(['ODE hidden-mode beta %.3f lambda %.9f, residual %.3e; ' ...
    'final modal amplitudes %.4g %.4g %.4g\n'],beta,real(lambda), ...
    eigenResidual,modalAmplitude(end,1),modalAmplitude(end,2), ...
    modalAmplitude(end,3));
end

function state=local_refine(phi,state,p)
for iteration=1:10000
    response=phi(state);
    residual=norm(response-state)/max(norm(state),eps);
    if residual<2e-13; break; end
    state=(1-p)*state+p*response;
end
fprintf('Refined hidden-mode fixed point in %d steps, residual %.3e.\n', ...
    iteration,residual);
end
