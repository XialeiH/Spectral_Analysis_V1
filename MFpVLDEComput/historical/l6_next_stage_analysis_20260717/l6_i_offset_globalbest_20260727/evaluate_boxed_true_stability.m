function evaluate_boxed_true_stability()
% Compare true moved-state stability at the boxed FPP and blue HC-minimum points.

setupFile = getenv('BOXED_SETUP_FILE');
orangeFile = getenv('BOXED_ORANGE_FILE');
blueFile = getenv('BOXED_BLUE_FILE');
outputRoot = getenv('BOXED_OUTPUT_ROOT');
if ~isfile(setupFile) || ~isfile(orangeFile) || ~isfile(blueFile) || isempty(outputRoot)
    error('BoxedStability:Environment','Required input paths are missing.');
end
if ~exist(outputRoot,'dir'); mkdir(outputRoot); end

loaded = load(setupFile,'setup');
setup = loaded.setup;
context = setup.Context;
orange = load(orangeFile,'fixed');
blue = load(blueFile,'fixedPoint');

beta6 = 0.034;
betaI = [0.0460053091161357;0.0176002942418064];
labels = ["orange FPP coordinate in true model";"blue HC-minimum coordinate"];
states = {orange.fixed.State;blue.fixedPoint};
baselineMaxReal = local_leading(setup.Pathway.JBaseline);
trueMaxReal = nan(2,1);
relativeChange = nan(2,1);
jvpRelativeError = nan(2,1);
for index = 1:2
    gain6 = 1+beta6;
    gainI = 1+betaI(index);
    state = states{index}(:);
    J = real_tuning_true_jacobian(state,context,gain6,gainI);
    trueMaxReal(index) = local_leading(J);
    relativeChange(index) = (trueMaxReal(index)-baselineMaxReal)/abs(baselineMaxReal);
    phi = @(x) l6ns_phi(x,1-gain6,context,[1 1],gainI,'true');
    jvpRelativeError(index) = local_jvp_error(J,phi,state,7300+index);
end

result = table(labels,repmat(beta6,2,1),betaI,repmat(baselineMaxReal,2,1), ...
    trueMaxReal,trueMaxReal-repmat(baselineMaxReal,2,1),relativeChange, ...
    jvpRelativeError,'VariableNames',{'caseLabel','beta6','betaI', ...
    'baselineMaxReal','trueMovedMaxReal','absoluteChange','relativeChange', ...
    'analyticJvpRelativeError'});
writetable(result,fullfile(outputRoot,'11_Boxed_TrueMoved_MaxReal_Comparison.tsv'), ...
    'FileType','text','Delimiter','\t');
save(fullfile(outputRoot,'11_Boxed_TrueMoved_MaxReal_Comparison.mat'), ...
    'result','-v7.3');
disp(result);
end

function maximumReal = local_leading(matrix)
options = struct('tol',1e-9,'maxit',1800,'p',80,'isreal',true,'disp',0);
try
    values = eigs(matrix,8,'largestreal',options);
catch
    values = eigs(matrix,8,'lr',options);
end
maximumReal = max(real(values));
end

function errorValue = local_jvp_error(J,phi,state,seed)
rng(seed,'twister');
direction = randn(size(state));
direction = direction/norm(direction);
step = 2e-6*max(1,norm(state)/sqrt(numel(state)));
numeric = (phi(state+step*direction)-phi(state-step*direction))/(2*step);
analytic = J*direction;
errorValue = norm(analytic-numeric)/max(norm(numeric),eps);
end
