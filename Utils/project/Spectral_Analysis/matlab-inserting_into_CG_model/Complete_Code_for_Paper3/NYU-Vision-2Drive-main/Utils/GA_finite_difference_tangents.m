function Dr = GA_finite_difference_tangents(responseMat, coordVals)
% Central finite-difference tangents for responses sampled along one coordinate.
% responseMat is nResponse-by-nCoord.
coordVals = coordVals(:)';
nCoord = numel(coordVals);
Dr = zeros(size(responseMat));

for k = 1:nCoord
    if k == 1
        h = coordVals(2) - coordVals(1);
        Dr(:,k) = (responseMat(:,2) - responseMat(:,1)) ./ h;
    elseif k == nCoord
        h = coordVals(end) - coordVals(end-1);
        Dr(:,k) = (responseMat(:,end) - responseMat(:,end-1)) ./ h;
    else
        h = coordVals(k+1) - coordVals(k-1);
        Dr(:,k) = (responseMat(:,k+1) - responseMat(:,k-1)) ./ h;
    end
end
end
