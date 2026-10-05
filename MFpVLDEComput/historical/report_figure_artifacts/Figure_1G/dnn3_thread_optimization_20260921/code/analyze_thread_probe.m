function analyze_thread_probe(outputRoot)
threads = [1 2 4 8 16];
seconds = zeros(10,5);
maxAbs = zeros(10,5); meanAbs = zeros(10,5);
sameIterations = false(10,5); residualDiff = zeros(10,5);
for repeat = 1:10
    ref = load(fullfile(outputRoot,'probe',sprintf( ...
        'dnn_03HC_threads_16_repeat_%03d.mat',repeat)),'state','T');
    reference = [ref.state.S;ref.state.C;ref.state.I];
    for k = 1:5
        d = load(fullfile(outputRoot,'probe',sprintf( ...
            'dnn_03HC_threads_%02d_repeat_%03d.mat',threads(k),repeat)),'state','T');
        assert(d.T.Converged && d.T.FinalResidual<5e-3 && d.T.Seconds>0);
        assert(d.T.NNThreads==threads(k) && d.T.AllocatedCPUs==16);
        current = [d.state.S;d.state.C;d.state.I];
        difference = abs(current-reference);
        assert(all(isfinite(difference)));
        seconds(repeat,k) = d.T.Seconds;
        maxAbs(repeat,k) = max(difference);
        meanAbs(repeat,k) = mean(difference);
        sameIterations(repeat,k) = d.T.Iterations==ref.T.Iterations;
        residualDiff(repeat,k) = abs(d.T.FinalResidual-ref.T.FinalResidual);
    end
end
eligible = all(maxAbs<=1e-3 & meanAbs<=1e-4 & sameIterations & residualDiff<=1e-6,1);
means = mean(seconds,1);
rank = means; rank(~eligible) = Inf;
[bestTime,k] = min(rank);
assert(isfinite(bestTime),'No numerically validated configuration.');
selectedThreads = threads(k);
summary = table(threads',means',std(seconds,0,1)', ...
    std(seconds,0,1)'/sqrt(10),max(maxAbs,[],1)',max(meanAbs,[],1)',eligible', ...
    'VariableNames',{'Threads','MeanSeconds','SD','SEM','MaxAbsDifference', ...
    'MaxMeanAbsDifference','Eligible'});
writetable(summary,fullfile(outputRoot,'probe_summary.csv'));
save(fullfile(outputRoot,'probe_summary.mat'),'seconds','summary','selectedThreads', ...
    'maxAbs','meanAbs','sameIterations','residualDiff');
fid=fopen(fullfile(outputRoot,'selected_threads.txt'),'w');
assert(fid>=0); fprintf(fid,'%d\n',selectedThreads); fclose(fid);
disp(summary);
fprintf('Selected %d threads. Probe mean %.9f s; fresh100-trial validation follows.\n', ...
    selectedThreads,bestTime);
end
