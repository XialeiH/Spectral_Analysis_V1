function run_narrow_instability_validation()
% Validate the narrow L6 unstable mode hidden from symmetric iteration.

paths = mechanism_initialize();
loadedSetup = load(paths.SetupFile,'setup');
context = loadedSetup.setup.Context;
loaded = load(fullfile(paths.OutputRoot,'branch_continuation.mat'), ...
    'l6ForwardStates');
betaGrid = 0:0.01:0.35;
beta = 0.17;
index = find(abs(betaGrid-beta)<1e-12,1);
fixed = loaded.l6ForwardStates(:,index);
previous = loaded.l6ForwardStates(:,index-1);
gain6 = 1+beta;
phi = @(x)mechanism_phi_variant(x,context,gain6,1,'extended');
J = real_tuning_true_jacobian(fixed,context,gain6,1);
[lambda,right,left,eigenResidual] = local_leading_pair(J);

stateScale = max(1,norm(fixed)/sqrt(numel(fixed)));
step = 2e-6*stateScale;
direction = real(right)/norm(real(right));
numeric = (phi(fixed+step*direction)-phi(fixed-step*direction))/(2*step);
jvpError = norm(J*direction-numeric)/max(norm(numeric),eps);
left = left/(left'*right);
continuationDisplacement = previous-fixed;
continuationProjection = abs(left'*continuationDisplacement)/ ...
    max(norm(left)*norm(continuationDisplacement),eps);

relaxation = context.RelaxationP;
relaxedMultiplier = 1-relaxation+relaxation*lambda;
epsilon = 1e-6*norm(fixed);
state = fixed+epsilon*direction;
maximumIterations = 1600;
relativeDistance = nan(maximumIterations+1,1);
modalAmplitude = nan(maximumIterations+1,1);
relativeDistance(1) = norm(state-fixed)/epsilon;
modalAmplitude(1) = abs(left'*(state-fixed))/epsilon;
for iteration = 1:maximumIterations
    state = (1-relaxation)*state+relaxation*phi(state);
    relativeDistance(iteration+1) = norm(state-fixed)/epsilon;
    modalAmplitude(iteration+1) = abs(left'*(state-fixed))/epsilon;
end
iterationGrid = (0:maximumIterations)';
linearPrediction = abs(relaxedMultiplier).^iterationGrid*modalAmplitude(1);

validation = table(beta,real(lambda),imag(lambda),abs(relaxedMultiplier), ...
    eigenResidual,jvpError,continuationProjection,relativeDistance(end), ...
    modalAmplitude(end), ...
    'VariableNames',{'beta','leadingReal','leadingImag', ...
    'relaxedIterationMultiplierMagnitude','eigenpairResidual', ...
    'finiteDifferenceJvpError','previousBranchStepProjectionOnUnstableLeftMode', ...
    'finalRelativeStateDistance','finalRelativeModalAmplitude'});
writetable(validation,fullfile(paths.OutputRoot, ...
    'narrow_instability_validation.tsv'),'FileType','text','Delimiter','\t');
save(fullfile(paths.OutputRoot,'narrow_instability_validation.mat'), ...
    'validation','iterationGrid','relativeDistance','modalAmplitude', ...
    'linearPrediction','lambda','right','left','-v7.3');

figure('Color','w','Position',[100 100 780 520]);
semilogy(iterationGrid,modalAmplitude,'LineWidth',1.6); hold on
semilogy(iterationGrid,linearPrediction,'--','LineWidth',1.5);
semilogy(iterationGrid,relativeDistance,':','LineWidth',1.5);
grid on; xlabel('relaxed iteration step');
ylabel('growth relative to initial perturbation');
title(sprintf(['L6 beta %.3f: unstable mode hidden from symmetric continuation, ' ...
    'lambda=%.6f'],beta,real(lambda)));
legend({'measured left-mode amplitude','linear prediction','total state distance'}, ...
    'Location','northwest');
exportgraphics(gcf,fullfile(paths.FigureRoot, ...
    '04_hidden_narrow_instability_perturbation_growth.pdf'),'ContentType','vector');
fprintf(['beta %.3f lambda %.9f%+.3gi, |mu_relaxed| %.9f, eig residual %.3e, ' ...
    'JVP %.3e, continuation projection %.3e, modal growth %.3g\n'], ...
    beta,real(lambda),imag(lambda),abs(relaxedMultiplier),eigenResidual, ...
    jvpError,continuationProjection,modalAmplitude(end));
end

function [lambda,right,left,residual] = local_leading_pair(J)
options = struct('tol',1e-11,'maxit',2600,'p',120,'isreal',true,'disp',0);
[rightVectors,rightValues] = eigs(J,12,'largestreal',options);
rightValues = diag(rightValues);
[~,rightIndex] = max(real(rightValues));
lambda = rightValues(rightIndex);
right = rightVectors(:,rightIndex);
residual = norm(J*right-lambda*right)/max(norm(J*right),eps);
[leftVectors,leftValues] = eigs(J',16,'largestreal',options);
leftValues = diag(leftValues);
[~,leftIndex] = min(abs(leftValues-conj(lambda)));
left = leftVectors(:,leftIndex);
end
