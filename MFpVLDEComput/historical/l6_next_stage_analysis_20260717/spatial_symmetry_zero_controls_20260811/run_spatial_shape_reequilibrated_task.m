function summary = run_spatial_shape_reequilibrated_task(taskIndex,setupFile,outputRoot)
% Re-equilibrate the h96 model after replacing the L4 spatial footprint.

if nargin<1 || isempty(taskIndex)
    taskIndex = str2double(getenv('SLURM_ARRAY_TASK_ID'));
end
if nargin<2 || isempty(setupFile)
    setupFile = getenv('SPATIAL_CONTROL_SETUP');
end
if nargin<3 || isempty(outputRoot)
    outputRoot = getenv('SPATIAL_CONTROL_OUTPUT');
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

controls = table(["left_half";"triangle"], ...
    ["L4_left_half_reequilibrated";"L4_triangle_reequilibrated"], ...
    'VariableNames',{'shape','control'});
if taskIndex<1 || taskIndex>height(controls)
    error('SpatialShape:Task','Task index must be 1 or 2.');
end

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
control = controls(taskIndex,:);
[operators,shapeDetails] = l6ns_spatial_shape_operators( ...
    context,control.shape,"baseline");

n = prod(context.MapSize);
initialVector = context.FixedPoint(:);
state = local_unpack(initialVector,n);
relaxation = 0.33;
blockEpochs = 250;
maximumEpochs = 3000;
tolerance = 1e-9;
residualHistory = nan(maximumEpochs/blockEpochs,1);
nanFlag = false;
elapsed = tic;

for block = 1:numel(residualHistory)
    [history,~,~,~,~,~,nanFlag] = local_iterate( ...
        state,relaxation,blockEpochs,context,operators);
    if nanFlag
        error('SpatialShape:IterationNaN','NaN during %s iteration.',control.control);
    end
    state = history{end};
    [oneStep,~,~,~,~,~,nanFlag] = local_iterate( ...
        state,1,1,context,operators);
    responseVector = local_pack(oneStep{end});
    stateVector = local_pack(state);
    residualHistory(block) = norm(responseVector-stateVector)/ ...
        max(norm(stateVector),eps);
    fprintf('%s epoch %d residual %.6e\n',control.control, ...
        block*blockEpochs,residualHistory(block));
    if residualHistory(block)<=tolerance
        break
    end
end
iterationTime = toc(elapsed);
iterations = block*blockEpochs;
residualHistory = residualHistory(1:block);
fixedPoint = local_pack(state);
converged = residualHistory(end)<=tolerance;
relativeBaselineDistance = norm(fixedPoint-initialVector)/ ...
    max(norm(initialVector),eps);

fprintf('%s: constructing Jacobian at its new equilibrium.\n',control.control);
[J,jacobianDetails] = l6ns_spatial_shape_control_jacobian( ...
    setup,control.shape,"baseline",fixedPoint,true);
eigenvalues = eig(full(J),'vector');
singularValues = svd(full(J));
e = 1:2*n;
i = 2*n+(1:n);
A = sparse(J(e,e)); B = sparse(J(e,i));
C = sparse(J(i,e)); D = sparse(J(i,i));
H = sparse(D-C*(A\B));
schurEigenvalues = eig(full(H),'vector');

summary = table(taskIndex,control.control,control.shape,iterations, ...
    converged,residualHistory(end),iterationTime,relativeBaselineDistance, ...
    max(real(eigenvalues)),min(real(eigenvalues)), ...
    sum(abs(eigenvalues)<=0.01),sum(abs(eigenvalues)<=0.025), ...
    sum(abs(eigenvalues)<=0.05),sum(singularValues<=0.05), ...
    sum(abs(schurEigenvalues)<=0.05), ...
    shapeDetails.L4RowSumRelativeError,shapeDetails.L4RetainedMassFraction, ...
    jacobianDetails.L4TranslationError,jacobianDetails.L4Rotation90Error, ...
    'VariableNames',{'taskIndex','control','l4Shape','iterations', ...
    'converged','fixedPointResidual','iterationSeconds', ...
    'relativeBaselineDistance','maximumRealEigenvalue','minimumRealEigenvalue', ...
    'countAbsEigLE0p01','countAbsEigLE0p025','countAbsEigLE0p05', ...
    'countSingularLE0p05','countSchurEigLE0p05', ...
    'l4RowSumRelativeError','l4RetainedMassFraction', ...
    'l4TranslationError','l4Rotation90Error'});

tag = char(control.control);
writetable(summary,fullfile(outputRoot,[tag '.tsv']), ...
    'FileType','text','Delimiter','\t');
equilibrium = state;
save(fullfile(outputRoot,[tag '.mat']),'summary','equilibrium','fixedPoint', ...
    'residualHistory','eigenvalues','singularValues','schurEigenvalues', ...
    'shapeDetails','-v7.3');
fprintf('%s: residual %.3e, eig/sigma/Schur <=.05 %d/%d/%d.\n',tag, ...
    summary.fixedPointResidual,summary.countAbsEigLE0p05, ...
    summary.countSingularLE0p05,summary.countSchurEigLE0p05);
end

function [history,a,b,c,d,e,nanFlag] = local_iterate( ...
        state,relaxation,epochs,context,operators)
placeholders = {[]};
emptyLibrary = {struct('S',{{}},'C',{{}},'I',{{}})};
[history,a,b,c,d,e,nanFlag] = ...
    LDEIteration_135FuncMain_CombDom_RealLGNL6_MLP6D_prefAngle( ...
    context.PixLGNCtgr,context.L6Kernel,state,relaxation, ...
    context.L6Parameters,epochs, ...
    operators.C_SS,operators.C_CS,operators.C_IS, ...
    operators.C_SC,operators.C_CC,operators.C_IC, ...
    operators.C_SI,operators.C_CI,operators.C_II, ...
    context.L4SEp,context.L4SIp,context.L4CEp,context.L4CIp, ...
    context.L4IEp,context.L4IIp,placeholders,placeholders,emptyLibrary, ...
    context.ContrastUse,context.OrientationUse,4,10,10, ...
    context.Isaturation,'xn',context.EKpUse,context.IKpUse);
end

function state = local_unpack(vector,n)
state = struct('S',vector(1:n),'C',vector(n+(1:n)),'I',vector(2*n+(1:n)));
end

function vector = local_pack(state)
vector = [state.S(:);state.C(:);state.I(:)];
end
