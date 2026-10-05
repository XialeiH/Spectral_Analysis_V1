function tangents = GA_resolvent_tangents(A, dPhi)
% Solve (I-A) dr = dPhi without forming inv(I-A).
n = size(A, 1);
tangents = (speye(n) - A) \ dPhi;
end
