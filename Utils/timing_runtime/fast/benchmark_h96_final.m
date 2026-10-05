function benchmark_h96_final()
% Final repeated benchmark for original, strict-fast, and float32-NN modes.
fastDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fastDir);
addpath(fullfile(rootDir, 'local_h96_manual_runtime'), '-begin');
addpath(fullfile(rootDir, 'Utils'), '-end');
addpath(fastDir, '-begin');
data = load(fullfile(fastDir, 'benchmark_input_c100_a0.mat'));
modelPath = fullfile(fastDir, 'h96_models_double.mat');

timer = tic;
baseDouble = h96_nn_prepare_base(data.ctx, 'double', modelPath);
baseDoubleFirstSeconds = toc(timer);
timer = tic;
strict = h96_nn_set_condition(baseDouble, data.Contr, data.angle);
strictConditionFirstSeconds = toc(timer);
timer = tic;
baseSingle = h96_nn_prepare_base(data.ctx, 'single', modelPath);
baseSingleFirstSeconds = toc(timer);
timer = tic;
rapid = h96_nn_set_condition(baseSingle, data.Contr, data.angle);
rapidConditionFirstSeconds = toc(timer);

% Warm every path before measuring interleaved repeats.
reference = run_reference(data);
strictState = h96_nn_fixedpoint_compact(data.IniTest, strict, 50);
rapidState = h96_nn_fixedpoint_compact(data.IniTest, rapid, 50);

repeatCount = 10;
referenceTimes = zeros(repeatCount, 1);
strictTimes = zeros(repeatCount, 1);
rapidTimes = zeros(repeatCount, 1);
for repeat = 1:repeatCount
    timer = tic;
    reference = run_reference(data);
    referenceTimes(repeat) = toc(timer);
    timer = tic;
    strictState = h96_nn_fixedpoint_compact(data.IniTest, strict, 50);
    strictTimes(repeat) = toc(timer);
    timer = tic;
    rapidState = h96_nn_fixedpoint_compact(data.IniTest, rapid, 50);
    rapidTimes(repeat) = toc(timer);
end

referenceVector = state_vector(reference);
strictVector = state_vector(strictState);
rapidVector = state_vector(rapidState);
result = struct();
result.Condition = struct('Contrast', data.Contr, 'Angle', data.angle, 'Iterations', 50);
result.ReferenceTimes = referenceTimes;
result.StrictTimes = strictTimes;
result.RapidTimes = rapidTimes;
result.ReferenceMedianSeconds = median(referenceTimes);
result.StrictMedianSeconds = median(strictTimes);
result.RapidMedianSeconds = median(rapidTimes);
result.StrictSpeedup = result.ReferenceMedianSeconds/result.StrictMedianSeconds;
result.RapidSpeedup = result.ReferenceMedianSeconds/result.RapidMedianSeconds;
result.StrictMaxAbsDifference = max(abs(strictVector-referenceVector));
result.StrictRelativeDifference = norm(strictVector-referenceVector)/norm(referenceVector);
result.RapidMaxAbsDifference = max(abs(rapidVector-referenceVector));
result.RapidRelativeDifference = norm(rapidVector-referenceVector)/norm(referenceVector);
result.BaseDoubleFirstSeconds = baseDoubleFirstSeconds;
result.BaseSingleFirstSeconds = baseSingleFirstSeconds;
result.StrictConditionFirstSeconds = strictConditionFirstSeconds;
result.RapidConditionFirstSeconds = rapidConditionFirstSeconds;
save(fullfile(fastDir, 'benchmark_h96_final.mat'), 'result', '-v7.3');

fprintf('\nH96_FINAL_50_ITERATION_BENCHMARK\n');
fprintf('reference_median_seconds=%.9f\n', result.ReferenceMedianSeconds);
fprintf('strict_median_seconds=%.9f\n', result.StrictMedianSeconds);
fprintf('strict_speedup=%.6f\n', result.StrictSpeedup);
fprintf('strict_max_abs_difference=%.12g\n', result.StrictMaxAbsDifference);
fprintf('strict_relative_difference=%.12g\n', result.StrictRelativeDifference);
fprintf('rapid_median_seconds=%.9f\n', result.RapidMedianSeconds);
fprintf('rapid_speedup=%.6f\n', result.RapidSpeedup);
fprintf('rapid_max_abs_difference=%.12g\n', result.RapidMaxAbsDifference);
fprintf('rapid_relative_difference=%.12g\n', result.RapidRelativeDifference);
fprintf('base_double_first_seconds=%.9f\n', baseDoubleFirstSeconds);
fprintf('base_single_first_seconds=%.9f\n', baseSingleFirstSeconds);
fprintf('strict_condition_first_seconds=%.9f\n', strictConditionFirstSeconds);
fprintf('rapid_condition_first_seconds=%.9f\n', rapidConditionFirstSeconds);
end

function state = run_reference(data)
ctx = data.ctx;
history = LDEIteration_135FuncMain_CombDom_RealLGNL6_MLP6D_prefAngle( ...
    ctx.PixLGNCtgr, ctx.L6Kernel, data.IniTest, ctx.p, ctx.L6pars{ctx.L6parId}, 50, ...
    ctx.C_SS_meanU, ctx.C_CS_meanU, ctx.C_IS_mean, ...
    ctx.C_SC_meanU, ctx.C_CC_meanU, ctx.C_IC_mean, ...
    ctx.C_SI_mean, ctx.C_CI_mean, ctx.C_II_mean, ...
    ctx.L4SEp, ctx.L4SIp, ctx.L4CEp, ctx.L4CIp, ctx.L4IEp, ctx.L4IIp, ...
    data.L4EmeshXAll, data.L4ImeshYAll, data.LDEFrfuncAll, ...
    data.Contruse, data.Orientationuse, ctx.N_HCOutY, ctx.NPixX, ctx.NPixY, ...
    ctx.Isaturation, 'xn', ctx.EKpUse, ctx.IKpUse);
state = history{end};
end

function vector = state_vector(state)
vector = [state.S; state.C; state.I];
end
