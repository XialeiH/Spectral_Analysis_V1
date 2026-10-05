function prepare_dense_finite_time_trajectory_task(figureIndex,setupFile,sourceRoot,outputRoot)
% Recreate one ODE trajectory on the 0.5 ms tangent-analysis grid.

arguments
    figureIndex (1,1) double {mustBeInteger,mustBePositive}
    setupFile (1,:) char
    sourceRoot (1,:) char
    outputRoot (1,:) char
end

sourceFile = local_source_file(sourceRoot,figureIndex);
setupLoaded = load(setupFile,'setup');
setup = setupLoaded.setup;
sourceLoaded = load(sourceFile,'output');
source = sourceLoaded.output;

tauMs = setup.Config.TauMs;
l6Weight = source.L6Weight;
initialState = source.InitialState(:);
sampleTimesMs = source.Times(source.Times<=500);
tangentStepMs = 0.5;
maximumHorizonMs = 100;
trajectoryTimesMs = 0:tangentStepMs:(sampleTimesMs(end)+maximumHorizonMs);

context = setup.Context;
phi = @(state) l6ns_phi(state,l6Weight,context,[1 1],1,'fpp');
rhs = @(~,state) (phi(state)-state)/tauMs;
options = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);
solution = ode45(rhs,[trajectoryTimesMs(1) trajectoryTimesMs(end)], ...
    initialState,options);
states = deval(solution,trajectoryTimesMs);

fixed = context.FixedPoint(:);
n = prod(double(context.MapSize(:).'));
jBaseline = sparse(l6ns_state_jacobian(setup,'L6',fixed,0));
[baselineOutput,~,baselineInput] = svds(jBaseline,1,'largest');
if real(baselineOutput'*baselineInput)<0
    baselineOutput = -baselineOutput;
end
baselineSingularOutputEMap = reshape( ...
    (1-context.CWeight)*real(baselineOutput(1:n))+ ...
    context.CWeight*real(baselineOutput(n+(1:n))), ...
    double(context.MapSize(:).'));

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end
outputFile = fullfile(outputRoot,sprintf('figure_%02d_trajectory.mat',figureIndex));
save(outputFile,'states','trajectoryTimesMs','sampleTimesMs','l6Weight', ...
    'baselineSingularOutputEMap','-v7.3');
fprintf('Figure %d: saved %d trajectory states through %.1f ms.\n', ...
    figureIndex,numel(trajectoryTimesMs),trajectoryTimesMs(end));
end

function sourceFile = local_source_file(sourceRoot,figureIndex)
relativePaths = { ...
    'outputs/7.1_Perturbation_Transient_and_Return_Baseline_L6_data.mat', ...
    'outputs/7.2_Perturbation_Transient_and_Return_Enhanced_L6_1p15x_data.mat', ...
    'outputs/7.4_Perturbation_Transient_and_Return_Baseline_L6_LargePerturbation_HC0p30_data.mat', ...
    'outputs/7.5_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC1_data.mat', ...
    'outputs/7.6_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC5_data.mat', ...
    'fixed_point_row_figures/7.6_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p00_fixedpointrow_data.mat', ...
    'fixed_point_row_figures/7.7_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p01_fixedpointrow_data.mat', ...
    'fixed_point_row_figures/7.8_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p02_fixedpointrow_data.mat', ...
    'fixed_direction_boundary_figures/7.9_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p03_data.mat', ...
    'fixed_direction_boundary_figures/7.10_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p10_data.mat', ...
    'fixed_direction_boundary_figures/7.11_Perturbation_Transient_and_Return_Baseline_L6_FixedDirection_HC6p50_data.mat', ...
    'outputs/7.7_Perturbation_Transient_and_Return_Baseline_L6_RandomPositive_HC10_data.mat'};
if figureIndex>numel(relativePaths)
    error('Figure index %d is out of range.',figureIndex);
end
sourceFile = fullfile(sourceRoot,relativePaths{figureIndex});
end
