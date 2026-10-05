function run_followup_branch_task()
% Compute FPP, direct-baseline, low-branch, and high-branch stability separately.

paths = followup_initialize();
loaded = load(paths.SetupFile,'setup');
context = loaded.setup.Context;
setup = loaded.setup;
baseline = context.FixedPoint(:);
specs = followup_specs('branch');
taskId = local_task_id();
spec = specs(taskId,:);
[gain6,gainI] = local_gains(spec.pathway,spec.beta);

low = followup_branch_state(paths,spec.pathway,spec.beta,"low",context);
high = followup_branch_state(paths,spec.pathway,spec.beta,"high",context);
if spec.pathway=="L6"
    jFpp = setup.Pathway.JRest+gain6*setup.Pathway.J6+setup.Pathway.JI;
else
    jFpp = setup.Pathway.JRest+setup.Pathway.J6+gainI*setup.Pathway.JI;
end
jDirect = real_tuning_true_jacobian(baseline,context,gain6,gainI);
jLow = real_tuning_true_jacobian(low,context,gain6,gainI);
jHigh = real_tuning_true_jacobian(high,context,gain6,gainI);
maxFpp = mechanism_max_real(jFpp);
maxDirect = mechanism_max_real(jDirect);
maxLow = mechanism_max_real(jLow);
maxHigh = mechanism_max_real(jHigh);

phi = @(x)mechanism_phi_variant(x,context,gain6,gainI,'extended');
lowResidual = norm(phi(low)-low)/max(norm(low),eps);
highResidual = norm(phi(high)-high)/max(norm(high),eps);
lowDomain = mechanism_domain_stats(low,context,gain6,gainI, ...
    spec.pathway,spec.beta,"low");
highDomain = mechanism_domain_stats(high,context,gain6,gainI, ...
    spec.pathway,spec.beta,"high");
[lowMean,lowMax] = local_population_stats(low);
[highMean,highMax] = local_population_stats(high);

row = table(taskId,spec.pathway,spec.beta,gain6,gainI,maxFpp,maxDirect, ...
    maxLow,maxHigh,lowResidual,highResidual, ...
    lowDomain.fractionOutsideTrainingBox,lowDomain.fractionFullExtension, ...
    highDomain.fractionOutsideTrainingBox,highDomain.fractionFullExtension, ...
    lowMean(1),lowMean(2),lowMean(3),lowMax(1),lowMax(2),lowMax(3), ...
    highMean(1),highMean(2),highMean(3),highMax(1),highMax(2),highMax(3), ...
    'VariableNames',{'taskId','pathway','beta','gain6','gainI', ...
    'maxRealFpp','maxRealDirectBaseline','maxRealLowBranch','maxRealHighBranch', ...
    'lowResidual','highResidual','lowFractionOutside','lowFractionFullExtension', ...
    'highFractionOutside','highFractionFullExtension', ...
    'lowMeanS','lowMeanC','lowMeanI','lowMaxS','lowMaxC','lowMaxI', ...
    'highMeanS','highMeanC','highMeanI','highMaxS','highMaxC','highMaxI'});
outDir = fullfile(paths.FollowupOutput,'branch_tasks');
if ~exist(outDir,'dir'); mkdir(outDir); end
writetable(row,fullfile(outDir,sprintf('branch_%03d.tsv',taskId)), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outDir,sprintf('branch_%03d.mat',taskId)), ...
    'row','low','high','-v7.3');
fprintf(['%s beta=%+.6f FPP/direct/low/high %.9f %.9f %.9f %.9f; ' ...
    'low/high max %.3f %.3f Hz\n'],spec.pathway,spec.beta,maxFpp,maxDirect, ...
    maxLow,maxHigh,max(low),max(high));
end

function taskId = local_task_id()
taskId = str2double(getenv('FOLLOWUP_TASK_ID'));
if ~isfinite(taskId); taskId = str2double(getenv('SLURM_ARRAY_TASK_ID')); end
if ~isfinite(taskId); taskId = 1; end
taskId = round(taskId);
end

function [gain6,gainI] = local_gains(pathway,beta)
gain6 = 1; gainI = 1;
if pathway=="L6"; gain6=1+beta; else; gainI=1+beta; end
end

function [means,maxima] = local_population_stats(state)
n = numel(state)/3;
values = [state(1:n),state(n+(1:n)),state(2*n+(1:n))];
means = mean(values,1);
maxima = max(values,[],1);
end
