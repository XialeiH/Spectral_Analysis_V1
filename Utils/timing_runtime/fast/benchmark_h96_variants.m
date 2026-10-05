function benchmark_h96_variants()
% Compare candidate h96 implementations from one reusable condition bundle.
fastDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fastDir);
addpath(fullfile(rootDir, 'local_h96_manual_runtime'), '-begin');
addpath(fullfile(rootDir, 'Utils'), '-end');
addpath(fastDir, '-begin');

data = load(fullfile(fastDir, 'benchmark_input_c100_a0.mat'));
fast = h96_nn_prepare_fast_config(data.ctx, data.Contr, data.angle, ...
    fullfile(fastDir, 'h96_models_double.mat'));
pageDouble = h96_nn_prepare_page_config(fast, 'double');
pageSingle = h96_nn_prepare_page_config(fast, 'single');
compactDouble = h96_nn_prepare_compact_config(pageDouble);
compactSingle = h96_nn_prepare_compact_config(pageSingle);

reference = run_reference(data);
referenceVector = state_vector(reference);

variants = { ...
    'pages_double', @() h96_nn_fixedpoint_pages(data.IniTest, pageDouble, 50); ...
    'compact_double', @() h96_nn_fixedpoint_compact(data.IniTest, compactDouble, 50); ...
    'compact_single_nn', @() h96_nn_fixedpoint_compact(data.IniTest, compactSingle, 50)};

repeatCount = 5;
results = struct([]);
fprintf('\nH96_VARIANT_BENCHMARK\n');
fprintf('%-23s %11s %10s %13s %13s\n', ...
    'variant', 'median_sec', 'speedup', 'max_abs_diff', 'relative_diff');
for variant = 1:size(variants, 1)
    name = variants{variant, 1};
    runner = variants{variant, 2};
    candidate = runner();
    times = zeros(repeatCount, 1);
    for repeat = 1:repeatCount
        timer = tic;
        candidate = runner();
        times(repeat) = toc(timer);
    end
    candidateVector = state_vector(candidate);
    results(variant).Name = name;
    results(variant).Times = times;
    results(variant).MedianSeconds = median(times);
    results(variant).MaxAbsDifference = max(abs(candidateVector - referenceVector));
    results(variant).RelativeDifference = norm(candidateVector - referenceVector) / norm(referenceVector);
    fprintf('%-23s %11.6f %10.3f %13.6g %13.6g\n', name, ...
        results(variant).MedianSeconds, ...
        results(1).MedianSeconds / results(variant).MedianSeconds, ...
        results(variant).MaxAbsDifference, results(variant).RelativeDifference);
end
save(fullfile(fastDir, 'benchmark_h96_variants.mat'), 'results', '-v7.3');
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
