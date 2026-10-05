function validate_figure1g_symmetry_50(runRoot)
% Check the accepted repeated-field accelerator over all 50 iterations.
addpath(fullfile(runRoot,'repo','Utils'),'-begin');
addpath(fullfile(runRoot,'fast'),'-begin');
addpath(fullfile(runRoot,'code'),'-begin');
threads = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isfinite(threads) && threads>=1
    maxNumCompThreads(threads);
end
modelFile = fullfile(runRoot,'fast','h96_models_double.mat');
for fieldHC = [4 40]
    B = load(fullfile(runRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)), ...
        'ctx','initialState');
    base = h96_nn_prepare_base(B.ctx,'single',modelFile);
    fastFull = h96_nn_set_condition_staticfirst(base,100,0);
    fastSymmetry = h96_nn_set_condition_symmetry(base,100,0);
    [fullState,fullInfo] = h96_nn_fixedpoint_staticfirst(B.initialState,fastFull,50);
    [symState,symInfo] = h96_nn_fixedpoint_symmetry(B.initialState,fastSymmetry,50);
    xFull = [fullState.S;fullState.C;fullState.I];
    xSym = [symState.S;symState.C;symState.I];
    fprintf(['FIELD=%d FULL50=%.6f SYMMETRY50=%.6f SPEEDUP=%.3f ', ...
        'RELERR=%.3e MAXERR=%.3e CHECKSUM_FULL=%.12g CHECKSUM_SYM=%.12g\n'], ...
        fieldHC,fullInfo.Seconds,symInfo.Seconds,fullInfo.Seconds/symInfo.Seconds, ...
        norm(xFull-xSym)/norm(xFull),max(abs(xFull-xSym)),sum(xFull),sum(xSym));
end
end
