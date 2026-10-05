function profile_h96_50iter()
% Profile the reference and isolated fast 50-step h96 fixed-point paths.

fastDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fastDir);
runtimeDir = fullfile(rootDir, 'local_h96_manual_runtime');
inputPath = fullfile(fastDir, 'benchmark_input_c100_a0.mat');
modelPath = fullfile(fastDir, 'h96_models_double.mat');

addpath(runtimeDir, '-begin');
addpath(fullfile(rootDir, 'Utils'), '-end');
addpath(fastDir, '-begin');

data = load(inputPath);
fast = h96_nn_prepare_fast_config(data.ctx, data.Contr, data.angle, modelPath);
pageDouble = h96_nn_prepare_page_config(fast, 'double');
pageSingle = h96_nn_prepare_page_config(fast, 'single');
compactDouble = h96_nn_prepare_compact_config(pageDouble);
compactSingle = h96_nn_prepare_compact_config(pageSingle);

profile clear;
profile on;
run_reference(data);
profile off;
referenceProfile = profile('info');
print_profile('REFERENCE', referenceProfile, 25);

profile clear;
profile on;
h96_nn_fixedpoint_fast(data.IniTest, fast, 50);
profile off;
fastProfile = profile('info');
print_profile('FAST', fastProfile, 25);

profile clear;
profile on;
h96_nn_fixedpoint_pages(data.IniTest, pageDouble, 50);
profile off;
pageDoubleProfile = profile('info');
print_profile('PAGES_DOUBLE', pageDoubleProfile, 25);

profile clear;
profile on;
h96_nn_fixedpoint_pages(data.IniTest, pageSingle, 50);
profile off;
pageSingleProfile = profile('info');
print_profile('PAGES_SINGLE', pageSingleProfile, 25);

profile clear;
profile on;
h96_nn_fixedpoint_compact(data.IniTest, compactDouble, 50);
profile off;
compactDoubleProfile = profile('info');
print_profile('COMPACT_DOUBLE', compactDoubleProfile, 25);

profile clear;
profile on;
h96_nn_fixedpoint_compact(data.IniTest, compactSingle, 50);
profile off;
compactSingleProfile = profile('info');
print_profile('COMPACT_SINGLE', compactSingleProfile, 25);

save(fullfile(fastDir, 'profile_h96_50iter.mat'), ...
    'referenceProfile', 'fastProfile', 'pageDoubleProfile', 'pageSingleProfile', ...
    'compactDoubleProfile', 'compactSingleProfile', '-v7.3');
end

function state = run_reference(data)
ctx = data.ctx;
history = LDEIteration_135FuncMain_CombDom_RealLGNL6_MLP6D_prefAngle( ...
    ctx.PixLGNCtgr, ctx.L6Kernel, data.IniTest, ctx.p, ...
    ctx.L6pars{ctx.L6parId}, 50, ...
    ctx.C_SS_meanU, ctx.C_CS_meanU, ctx.C_IS_mean, ...
    ctx.C_SC_meanU, ctx.C_CC_meanU, ctx.C_IC_mean, ...
    ctx.C_SI_mean, ctx.C_CI_mean, ctx.C_II_mean, ...
    ctx.L4SEp, ctx.L4SIp, ctx.L4CEp, ctx.L4CIp, ctx.L4IEp, ctx.L4IIp, ...
    data.L4EmeshXAll, data.L4ImeshYAll, data.LDEFrfuncAll, ...
    data.Contruse, data.Orientationuse, ...
    ctx.N_HCOutY, ctx.NPixX, ctx.NPixY, ctx.Isaturation, ...
    'xn', ctx.EKpUse, ctx.IKpUse);
state = history{end};
end

function print_profile(label, profileData, rowCount)
tableData = profileData.FunctionTable;
times = [tableData.TotalTime];
[~, order] = sort(times, 'descend');
fprintf('\n%s_PROFILE\n', label);
fprintf('%10s %10s  %s\n', 'seconds', 'calls', 'function');
for index = order(1:min(rowCount, numel(order)))
    item = tableData(index);
    fprintf('%10.6f %10d  %s\n', item.TotalTime, item.NumCalls, item.FunctionName);
end
end
