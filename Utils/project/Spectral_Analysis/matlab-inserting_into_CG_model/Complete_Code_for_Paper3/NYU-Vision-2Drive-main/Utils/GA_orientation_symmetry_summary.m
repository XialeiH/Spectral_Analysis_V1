function out = GA_orientation_symmetry_summary(responseMat, anglesDeg)
% Orientation/direction consistency summaries for section 7.
% responseMat is nResponse-by-nAngle.
anglesDeg = anglesDeg(:)';
nAngle = numel(anglesDeg);
out.AnglesDeg = anglesDeg;
out.OrientationPeriodError = NaN;
out.DirectionClosureError = NaN;
out.NeighborStepNorm = zeros(1, max(nAngle - 1, 0));

for k = 1:nAngle-1
    out.NeighborStepNorm(k) = norm(responseMat(:,k+1) - responseMat(:,k));
end

idx180 = find(abs(anglesDeg - (anglesDeg(1) + 180)) < 1e-8, 1);
if ~isempty(idx180)
    out.OrientationPeriodError = norm(responseMat(:,idx180) - responseMat(:,1));
end

idx360 = find(abs(anglesDeg - (anglesDeg(1) + 360)) < 1e-8, 1);
if ~isempty(idx360)
    out.DirectionClosureError = norm(responseMat(:,idx360) - responseMat(:,1));
end
end
