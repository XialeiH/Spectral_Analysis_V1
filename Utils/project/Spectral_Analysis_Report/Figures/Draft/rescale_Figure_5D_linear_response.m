function rescale_Figure_5D_linear_response
% Exact amplitude conversion for the stored homogeneous linear ODE solutions.
root = [repro_paths('project') '/Spectral_Analysis/matlab-inserting_into_CG_model/report_figure_artifacts/6_Pathway_Compensation/Figure_5D/fixed_input_dynamic_difference'];
baselineFile = fullfile(root,'baseline.mat');
b = load(baselineFile);
scale = 1/b.initialHC;
b.initialDelta = scale*b.initialDelta;
b.baselineDelta = scale*b.baselineDelta;
b.baselineEnergy = scale^2*b.baselineEnergy;
b.initialHC = 1;
b.amplitudeMethod = 'Exact linear homogeneity rescaling to initial HC = 1 sp/s';
assert(abs(hc(b.initialDelta)-1)<1e-12);
save(baselineFile,'-struct','b');
metrics = readtable(fullfile(root,'dynamic_difference.tsv'), ...
    'FileType','text','Delimiter','\t');
for k = 1:height(metrics)
    file = fullfile(root,sprintf('point_%02d.mat',metrics.taskId(k)));
    p = load(file);
    before = hc(p.delta);
    [~,oldPeak] = max(before);
    scale = 1/p.initialHC;
    p.delta = scale*p.delta;
    p.initialDelta = scale*p.initialDelta;
    p.initialHC = 1;
    p.amplitudeMethod = b.amplitudeMethod;
    after = hc(p.delta);
    [~,newPeak] = max(after);
    assert(oldPeak==newPeak,'Amplitude scaling unexpectedly changed peak time.');
    assert(abs(after(1)-1)<1e-12);
    assert(isequal(p.initialDelta,b.initialDelta));
    discrepancy = sqrt(trapz(p.timesMs,hc(p.delta-b.baselineDelta).^2)/b.baselineEnergy);
    assert(abs(discrepancy-metrics.dynamicDifference(k))<1e-10);
    save(file,'-struct','p');
    fprintf('Task %02d: initial HC %.12g sp/s; peak HC %.8g sp/s at %.1f ms\n', ...
        metrics.taskId(k),after(1),max(after),p.timesMs(newPeak));
end
end

function value = hc(delta)
n = size(delta,1)/3;
e = 0.6923*delta(1:n,:)+0.3077*delta(n+(1:n),:);
i = delta(2*n+(1:n),:);
value = sqrt(mean(0.8*e.^2+0.2*i.^2,1));
end
