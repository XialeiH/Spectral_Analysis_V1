function run_figure1g_trial(runRoot, trialId)
% Run one paired CG/DNN timing trial and save timing scalars only.
trialId = str2double(string(trialId));
iterations = 100;
repoRoot = fullfile(runRoot, 'repo');
addpath(fullfile(repoRoot, 'Utils'), '-begin');
addpath(fullfile(runRoot, 'code'), '-begin');
addpath(fullfile(runRoot, 'fast'), '-begin');
maxNumCompThreads(1);

manifest = readtable(fullfile(runRoot, 'manifest.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t', 'TextType', 'string');
row = manifest(manifest.Trial == trialId, :);
assert(height(row) == 1, 'Expected exactly one manifest row for trial %d.', trialId);

bundle = load(fullfile(runRoot, 'fast', 'benchmark_input_c100_a0.mat'), 'ctx', 'IniTest');
ctx = bundle.ctx;
initialState = bundle.IniTest;
angle = row.Angle;
contrast = row.Contrast;
libraryFile = fullfile(repoRoot, 'Data', 'Paper2_NetworkTuning', 'Fig1V4', ...
    'Paper3PlotingData', row.LibraryFile);
assert(isfile(libraryFile), 'Missing library file: %s', libraryFile);

ctgrOrderReadout = paper3_readout_order(angle);
modelFile = fullfile(runRoot, 'fast', 'h96_models_double.mat');

% Alternate execution order to avoid a systematic first-method bias.
if mod(trialId, 2) == 1
    cg = cg_fixedpoint_timing_only(libraryFile, initialState, ctx, ...
        ctgrOrderReadout, iterations);
    dnn = run_dnn(initialState, ctx, modelFile, contrast, angle, iterations);
else
    dnn = run_dnn(initialState, ctx, modelFile, contrast, angle, iterations);
    cg = cg_fixedpoint_timing_only(libraryFile, initialState, ctx, ...
        ctgrOrderReadout, iterations);
end

assert(all(isfinite([cg.TotalSeconds, cg.SecondsPerIteration, ...
    dnn.TotalSeconds, dnn.SecondsPerIteration, cg.FinalChecksum, dnn.FinalChecksum])));

result = table(trialId, angle, contrast, iterations, ...
    cg.SecondsPerIteration, dnn.SecondsPerIteration, ...
    cg.TotalSeconds, dnn.TotalSeconds, cg.SetupSeconds, dnn.SetupSeconds, ...
    cg.FinalChecksum, dnn.FinalChecksum, ...
    'VariableNames', {'Trial','Angle','Contrast','Iterations', ...
    'CGSecondsPerIteration','DNNSecondsPerIteration', ...
    'CGTotalSeconds','DNNTotalSeconds','CGSetupSeconds','DNNSetupSeconds', ...
    'CGChecksum','DNNChecksum'});
outputFile = fullfile(runRoot, 'results', sprintf('trial_%03d.tsv', trialId));
writetable(result, outputFile, 'FileType', 'text', 'Delimiter', '\t');
fprintf('FIGURE1G trial=%d angle=%.2f contrast=%g cg_total=%.6f dnn_total=%.6f\n', ...
    trialId, angle, contrast, cg.TotalSeconds, dnn.TotalSeconds);
end

function info = run_dnn(initialState, ctx, modelFile, contrast, angle, iterations)
timer = tic;
base = h96_nn_prepare_base(ctx, 'double', modelFile);
fast = h96_nn_set_condition(base, contrast, angle);
setupSeconds = toc(timer);
[state, loopInfo] = h96_nn_fixedpoint_compact(initialState, fast, iterations);
info.SetupSeconds = setupSeconds;
info.SecondsPerIteration = loopInfo.SecondsPerIteration;
info.TotalSeconds = toc(timer);
info.FinalChecksum = sum(state.S) + sum(state.C) + sum(state.I);
end

function order = paper3_readout_order(angle)
rotInd = floor(angle/45);
mirInd = mod(floor(angle/22.5), 2);
categoryOrder = [2 1; 3 4];
categoryOrder = rot90(categoryOrder, rotInd);
if mirInd
    categoryOrder = flip(categoryOrder, 1);
end
order = [categoryOrder(1,2); categoryOrder(1,1); ...
    categoryOrder(2,1); categoryOrder(2,2)];
end
