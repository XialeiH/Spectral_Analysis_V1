function audit = l6ei_first_order_slope(setupFile,outputRoot)
% Compute and finite-difference-check the FPP beta6-to-betaEI slope.

loaded=load(setupFile,'setup');
pathway=loaded.setup.Pathway;
jBaseline=pathway.JBaseline;
j6=pathway.J6;
jI=pathway.JI;
n=pathway.PopulationSize;
iColumns=2*n+(1:n);
jEI=sparse(size(jI,1),size(jI,2));
jEI(1:2*n,iColumns)=jI(1:2*n,iColumns);

[lambda,left,right]=local_leading_mode(jBaseline,loaded.setup.Config);
sensitivity6=real(left'*(j6*right));
sensitivityEI=real(left'*(jEI*right));
theorySlope=-sensitivity6/sensitivityEI;
step=1e-5;
finiteDifference6=local_central_sensitivity( ...
    jBaseline,j6,step,loaded.setup.Config);
finiteDifferenceEI=local_central_sensitivity( ...
    jBaseline,jEI,step,loaded.setup.Config);
finiteDifferenceSlope=-finiteDifference6/finiteDifferenceEI;

audit=table(real(lambda),imag(lambda),sensitivity6,sensitivityEI, ...
    finiteDifference6,finiteDifferenceEI,theorySlope, ...
    finiteDifferenceSlope, ...
    abs(finiteDifference6-sensitivity6)/max(abs(sensitivity6),eps), ...
    abs(finiteDifferenceEI-sensitivityEI)/max(abs(sensitivityEI),eps), ...
    norm(jEI,'fro'),norm(jI-jEI,'fro'),step, ...
    'VariableNames',{'baselineLeadingReal','baselineLeadingImag', ...
    'sensitivityL6','sensitivityItoE','finiteDifferenceL6', ...
    'finiteDifferenceItoE','slopeL6ItoE','finiteDifferenceSlope', ...
    'relativeErrorL6','relativeErrorItoE','ItoEJacobianFrobeniusNorm', ...
    'excludedItoIJacobianFrobeniusNorm','finiteDifferenceStep'});
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end
writetable(audit,fullfile(outputRoot,'theoretical_fpp_slope_audit.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputRoot,'theoretical_fpp_slope_audit.mat'), ...
    'audit','jEI','-v7.3');
fprintf(['L6-I-to-E FPP slope %.12g; finite-difference slope %.12g; ' ...
    'relative sensitivity errors %.3e and %.3e.\n'], ...
    theorySlope,finiteDifferenceSlope,audit.relativeErrorL6, ...
    audit.relativeErrorItoE);
end

function sensitivity=local_central_sensitivity(jBaseline,component,step,cfg)
alphaPlus=real(local_leading_mode(jBaseline+step*component,cfg));
alphaMinus=real(local_leading_mode(jBaseline-step*component,cfg));
sensitivity=(alphaPlus-alphaMinus)/(2*step);
end

function [lambda,left,right]=local_leading_mode(matrix,cfg)
count=6;
options=struct('tol',cfg.EigsTolerance,'maxit',cfg.EigsMaxIterations, ...
    'p',min(max(cfg.EigsSubspaceDimension,2*count+8),size(matrix,1)), ...
    'disp',0,'isreal',true);
try
    [rightPool,rightValues]=eigs(matrix,count,'largestreal',options);
catch
    [rightPool,rightValues]=eigs(matrix,count,'lr',options);
end
rightLambda=diag(rightValues);
[~,index]=max(real(rightLambda));
lambda=rightLambda(index);
right=rightPool(:,index);
right=right/norm(right);
try
    [leftPool,leftValues]=eigs(matrix',count,'largestreal',options);
catch
    [leftPool,leftValues]=eigs(matrix',count,'lr',options);
end
leftLambda=diag(leftValues);
[~,leftIndex]=min(abs(leftLambda-conj(lambda)));
left=leftPool(:,leftIndex);
overlap=left'*right;
if abs(overlap)<1e-12
    error('L6EI:Eigenvectors','Leading left/right overlap is singular.');
end
left=left/conj(overlap);
end
