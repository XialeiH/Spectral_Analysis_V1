function transient = l6ns_continuous_transient(jacobian, timeGrid, cfg, weightVector, tau)
% Top singular input/output pair of exp(t*(J-I)/tau) over a time grid.

n = size(jacobian, 1);
if nargin < 5 || isempty(tau)
    tau = 1;
end
if ~isscalar(tau) || ~isfinite(tau) || tau <= 0
    error('tau must be one positive finite scalar.');
end
generator = (jacobian - speye(n)) / tau;
if nargin<4 || isempty(weightVector)
    weightVector=ones(n,1);
    normType='unweighted';
else
    weightVector=weightVector(:);
    if numel(weightVector)~=n || any(~isfinite(weightVector)) || any(weightVector<=0)
        error('weightVector must contain one positive finite weight per state variable.');
    end
    normType='fixed_point_rms_scaled';
end
weightMatrix=spdiags(weightVector,0,n,n);
inverseWeightMatrix=spdiags(1./weightVector,0,n,n);
generator=weightMatrix*generator*inverseWeightMatrix;
rng(cfg.RandomSeed);
seed = randn(n,1);
seed = seed / norm(seed);

gain = nan(numel(timeGrid), 1);
inputVectors = cell(numel(timeGrid), 1);
outputVectors = cell(numel(timeGrid), 1);
for ti = 1:numel(timeGrid)
    time = timeGrid(ti);
    v = seed;
    for iteration = 1:cfg.TransientPowerIterations
        uRaw = l6ns_expmv_krylov(generator, time, v, 64, 1e-11);
        sigma = norm(uRaw);
        if sigma == 0; break; end
        u = uRaw / sigma;
        vRaw = l6ns_expmv_krylov(generator', time, u, 64, 1e-11);
        if norm(vRaw) == 0; break; end
        vNew = vRaw / norm(vRaw);
        if abs(vNew' * v) > 1 - 1e-8
            v = vNew;
            break
        end
        v = vNew;
    end
    uRaw = l6ns_expmv_krylov(generator, time, v, 72, 1e-11);
    gain(ti) = norm(uRaw);
    inputVectors{ti} = v;
    outputVectors{ti} = uRaw / max(gain(ti), eps);
    seed = v;
    fprintf('Transient time %.4g: gain %.8g.\n', time, gain(ti));
end

[peakGain, peakIndex] = max(gain);
peakTime = timeGrid(peakIndex);
svdsFlag = nan;
if isfield(cfg,'TransientRefinePeak') && cfg.TransientRefinePeak
    operator = @(x, flag) local_expm_multiply(generator,peakTime,x,flag);
    [uRefined,sRefined,vRefined,svdsFlag] = svds(operator,[n n],1,'largest', ...
        'Tolerance',cfg.TransientSvdsTolerance, ...
        'MaxIterations',cfg.TransientSvdsMaxIterations, ...
        'SubspaceDimension',cfg.TransientSvdsSubspaceDimension, ...
        'RightStartVector',inputVectors{peakIndex}, ...
        'FailureTreatment','keep');
    if svdsFlag~=0 || ~isfinite(sRefined(1,1))
        error('Peak-time svds refinement failed with flag %d.',svdsFlag);
    end
    peakGain = sRefined(1,1);
    gain(peakIndex) = peakGain;
    inputVectors{peakIndex} = vRefined(:,1);
    outputVectors{peakIndex} = uRefined(:,1);
end
transient = struct();
transient.TimeGrid = timeGrid(:);
transient.Gain = gain;
transient.PeakTime = peakTime;
transient.PeakGain = peakGain;
transient.OptimalInput = inverseWeightMatrix*inputVectors{peakIndex};
transient.PeakOutput = inverseWeightMatrix*outputVectors{peakIndex};
transient.WeightVector = weightVector;
transient.NormType = normType;
transient.Tau = tau;
transient.SingularSolverFlag = svdsFlag;

vScaled = inputVectors{peakIndex};
uScaled = outputVectors{peakIndex};
forward = l6ns_expmv_krylov(generator,transient.PeakTime,vScaled,72,1e-11);
adjoint = l6ns_expmv_krylov(generator',transient.PeakTime,uScaled,72,1e-11);
transient.ForwardSingularResidual = norm(forward-peakGain*uScaled) / ...
    max(norm(forward),eps);
transient.AdjointSingularResidual = norm(adjoint-peakGain*vScaled) / ...
    max(norm(adjoint),eps);
end

function y = local_expm_multiply(generator,time,x,flag)
if strcmp(flag,'notransp')
    y = l6ns_expmv_krylov(generator,time,x,72,1e-11);
else
    y = l6ns_expmv_krylov(generator',time,x,72,1e-11);
end
end
