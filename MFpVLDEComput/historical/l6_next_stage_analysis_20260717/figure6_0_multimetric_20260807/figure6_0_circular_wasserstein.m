function [circularDistance,linearDistance] = figure6_0_circular_wasserstein( ...
        baselineCurve,movedCurve,coarseAngles,fineStep)
% Compute normalized orientation-distribution W1 on a 180-degree circle.

if nargin<4; fineStep=0.1; end
uniqueMask = coarseAngles<180;
angles = coarseAngles(uniqueMask);
baseline = max(real(baselineCurve(uniqueMask)),0);
moved = max(real(movedCurve(uniqueMask)),0);
fineAngles = 0:fineStep:(180-fineStep);
baselineFine = interp1([angles(:);180],[baseline(:);baseline(1)], ...
    fineAngles,'pchip');
movedFine = interp1([angles(:);180],[moved(:);moved(1)], ...
    fineAngles,'pchip');
baselineFine = max(baselineFine,0);
movedFine = max(movedFine,0);
if sum(baselineFine)<=0 || sum(movedFine)<=0
    circularDistance=NaN;
    linearDistance=NaN;
    return
end
p = baselineFine/sum(baselineFine);
q = movedFine/sum(movedFine);
cumulativeDifference = cumsum(p-q);
linearDistance = fineStep*sum(abs(cumulativeDifference));
circularDistance = fineStep*sum(abs( ...
    cumulativeDifference-median(cumulativeDifference)));
end
