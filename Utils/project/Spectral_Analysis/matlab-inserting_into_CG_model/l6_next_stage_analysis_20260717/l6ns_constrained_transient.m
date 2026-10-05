function transient = l6ns_constrained_transient(jacobian,timeGridMs,cfg,operators)
% Compute C exp(Lt) B and its leading constrained singular channel.

n = size(jacobian,1);
generator = (jacobian-speye(n))/cfg.TauMs;
nTime = numel(timeGridMs);
nControl = size(operators.B,2);
nObservation = size(operators.C,1);
gain = nan(nTime,1);
kernels = nan(nObservation,nControl,nTime);
right = cell(nTime,1);
left = cell(nTime,1);
singularValues = cell(nTime,1);

for ti = 1:nTime
    response = zeros(n,nControl);
    for bi = 1:nControl
        response(:,bi) = l6ns_expmv_krylov( ...
            generator,timeGridMs(ti),operators.B(:,bi),72,1e-11);
    end
    h = full(operators.C*response);
    kernels(:,:,ti) = h;
    [u,s,v] = svd(h,'econ');
    gain(ti) = s(1,1);
    right{ti} = v(:,1);
    left{ti} = u(:,1);
    singularValues{ti} = diag(s);
end

[peakGain,peakIndex] = max(gain);
transient = struct();
transient.TimeGridMs = timeGridMs(:);
transient.Gain = gain;
transient.PeakGain = peakGain;
transient.PeakTimeMs = timeGridMs(peakIndex);
transient.OptimalControl = right{peakIndex};
transient.OptimalObservation = left{peakIndex};
transient.PeakSingularValues = singularValues{peakIndex};
positiveTime = find(timeGridMs>0);
if isempty(positiveTime)
    transient.PeakAfterZeroGain = peakGain;
    transient.PeakAfterZeroTimeMs = timeGridMs(peakIndex);
else
    [transient.PeakAfterZeroGain,relativeIndex] = max(gain(positiveTime));
    transient.PeakAfterZeroTimeMs = timeGridMs(positiveTime(relativeIndex));
end
transient.Kernels = kernels;
transient.ControlNames = operators.ControlNames;
transient.ObservationNames = operators.ObservationNames;
end
