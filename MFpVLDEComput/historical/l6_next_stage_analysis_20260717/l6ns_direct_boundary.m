function boundary = l6ns_direct_boundary(a, b, cfg)
% Locate max Re(lambda(A+(1-w)B)) = 1 by safeguarded bisection.

zeroW = 0;
zeroValue = l6ns_max_real(a + b, cfg);
if zeroValue < cfg.CriticalBoundary
    highW = zeroW;
    highValue = zeroValue;
    lowW = -0.2;
    lowValue = l6ns_max_real(a + (1 - lowW) * b, cfg);
    while lowValue < cfg.CriticalBoundary && lowW > -20
        highW = lowW;
        highValue = lowValue;
        lowW = 1.75 * lowW;
        lowValue = l6ns_max_real(a + (1 - lowW) * b, cfg);
    end
    if lowValue < cfg.CriticalBoundary
        error('Could not bracket the stability boundary down to w=%.6g.', lowW);
    end
elseif zeroValue > cfg.CriticalBoundary
    lowW = zeroW;
    lowValue = zeroValue;
    highW = 0.2;
    highValue = l6ns_max_real(a + (1 - highW) * b, cfg);
    while highValue > cfg.CriticalBoundary && highW < 20
        lowW = highW;
        lowValue = highValue;
        highW = 1.75 * highW;
        highValue = l6ns_max_real(a + (1 - highW) * b, cfg);
    end
    if highValue > cfg.CriticalBoundary
        error('Could not bracket the stability boundary up to w=%.6g.', highW);
    end
else
    [maxReal,lambda,rightVector] = l6ns_max_real(a + b, cfg);
    boundary = struct('W',zeroW,'Alpha',1,'MaxReal',maxReal, ...
        'Lambda',lambda,'RightVector',rightVector,'Bracket',[zeroW zeroW], ...
        'BracketValues',[zeroValue zeroValue],'Iterations',0);
    return
end

for iteration = 1:34
    middleW = 0.5 * (lowW + highW);
    middleValue = l6ns_max_real(a + (1 - middleW) * b, cfg);
    if middleValue >= cfg.CriticalBoundary
        lowW = middleW;
        lowValue = middleValue;
    else
        highW = middleW;
        highValue = middleValue;
    end
    if abs(highW - lowW) < 2e-8
        break
    end
end

boundary = struct();
boundary.W = 0.5 * (lowW + highW);
boundary.Alpha = 1 - boundary.W;
[boundary.MaxReal, boundary.Lambda, boundary.RightVector] = ...
    l6ns_max_real(a + boundary.Alpha * b, cfg);
boundary.Bracket = [lowW highW];
boundary.BracketValues = [lowValue highValue];
boundary.Iterations = iteration;
end
