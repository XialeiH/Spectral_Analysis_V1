function validate_h96_compact_conditions()
% Verify compact h96 solvers across low/high contrast and orientation inputs.
fastDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fastDir);
addpath(fullfile(rootDir, 'local_h96_manual_runtime'), '-begin');
addpath(fullfile(rootDir, 'Utils'), '-end');
addpath(fastDir, '-begin');
data = load(fullfile(fastDir, 'benchmark_input_c100_a0.mat'));
modelPath = fullfile(fastDir, 'h96_models_double.mat');

baseTimer = tic;
baseDouble = h96_nn_prepare_base(data.ctx, 'double', modelPath);
doubleBaseSeconds = toc(baseTimer);
baseTimer = tic;
baseSingle = h96_nn_prepare_base(data.ctx, 'single', modelPath);
singleBaseSeconds = toc(baseTimer);

conditions = [5 0; 42 7.5; 71 15; 100 22.5];
rows = struct([]);
fprintf('\nH96_COMPACT_MULTI_CONDITION_VALIDATION\n');
fprintf('%8s %8s %10s %10s %12s %12s %12s %12s\n', ...
    'contrast', 'angle', 'ref_sec', 'dbl_sec', 'dbl_max_abs', 'dbl_relative', ...
    'sgl_max_abs', 'sgl_relative');
for condition = 1:size(conditions, 1)
    contrast = conditions(condition, 1);
    angle = conditions(condition, 2);
    contrastVector = contrast*ones(baseDouble.n, 1);
    angleVector = angle*ones(size(contrastVector));

    timer = tic;
    reference = run_reference(data, contrastVector, angleVector);
    referenceSeconds = toc(timer);
    fastDouble = h96_nn_set_condition(baseDouble, contrast, angle);
    fastSingle = h96_nn_set_condition(baseSingle, contrast, angle);
    timer = tic;
    doubleState = h96_nn_fixedpoint_compact(data.IniTest, fastDouble, 50);
    doubleSeconds = toc(timer);
    timer = tic;
    singleState = h96_nn_fixedpoint_compact(data.IniTest, fastSingle, 50);
    singleSeconds = toc(timer);

    referenceVector = state_vector(reference);
    doubleVector = state_vector(doubleState);
    singleVector = state_vector(singleState);
    rows(condition).Contrast = contrast;
    rows(condition).Angle = angle;
    rows(condition).ReferenceSeconds = referenceSeconds;
    rows(condition).DoubleSeconds = doubleSeconds;
    rows(condition).SingleSeconds = singleSeconds;
    rows(condition).DoubleMaxAbs = max(abs(doubleVector-referenceVector));
    rows(condition).DoubleRelative = norm(doubleVector-referenceVector)/norm(referenceVector);
    rows(condition).SingleMaxAbs = max(abs(singleVector-referenceVector));
    rows(condition).SingleRelative = norm(singleVector-referenceVector)/norm(referenceVector);
    fprintf('%8g %8g %10.6f %10.6f %12.5g %12.5g %12.5g %12.5g\n', ...
        contrast, angle, referenceSeconds, doubleSeconds, rows(condition).DoubleMaxAbs, ...
        rows(condition).DoubleRelative, rows(condition).SingleMaxAbs, rows(condition).SingleRelative);
end
save(fullfile(fastDir, 'validation_h96_compact_conditions.mat'), ...
    'rows', 'doubleBaseSeconds', 'singleBaseSeconds', 'conditions', '-v7.3');
fprintf('double_base_prepare_seconds=%.6f\n', doubleBaseSeconds);
fprintf('single_base_prepare_seconds=%.6f\n', singleBaseSeconds);
end

function state = run_reference(data, contrastVector, angleVector)
ctx = data.ctx;
history = LDEIteration_135FuncMain_CombDom_RealLGNL6_MLP6D_prefAngle( ...
    ctx.PixLGNCtgr, ctx.L6Kernel, data.IniTest, ctx.p, ctx.L6pars{ctx.L6parId}, 50, ...
    ctx.C_SS_meanU, ctx.C_CS_meanU, ctx.C_IS_mean, ...
    ctx.C_SC_meanU, ctx.C_CC_meanU, ctx.C_IC_mean, ...
    ctx.C_SI_mean, ctx.C_CI_mean, ctx.C_II_mean, ...
    ctx.L4SEp, ctx.L4SIp, ctx.L4CEp, ctx.L4CIp, ctx.L4IEp, ctx.L4IIp, ...
    data.L4EmeshXAll, data.L4ImeshYAll, data.LDEFrfuncAll, ...
    contrastVector, angleVector, ctx.N_HCOutY, ctx.NPixX, ctx.NPixY, ...
    ctx.Isaturation, 'xn', ctx.EKpUse, ctx.IKpUse);
state = history{end};
end

function vector = state_vector(state)
vector = [state.S; state.C; state.I];
end
