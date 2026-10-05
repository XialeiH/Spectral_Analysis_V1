function dPhi = GA_dPhi_from_response_tangent(A, Dr)
% From (I-A)dr=dPhi, recover the feedforward task derivative dPhi.
dPhi = (speye(size(A, 1)) - A) * Dr;
end
