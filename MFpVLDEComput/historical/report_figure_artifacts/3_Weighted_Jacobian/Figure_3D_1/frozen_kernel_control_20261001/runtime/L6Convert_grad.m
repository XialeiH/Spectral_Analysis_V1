function dy = L6Convert_grad(x, L6pars)
% Derivative of the active L6Convert curve with respect to x.
% This keeps Jacobian construction aligned with all L6Convert modes,
% including nested/weighted parameter cells used by L6ParamUse.
if iscell(L6pars{end}) && numel(L6pars{end}) >= 3 && strcmp(L6pars{end}{1}, 'c1smooth')
    basePars = L6pars{end}{2};
    grid = L6pars{end}{3};
    yGrid = L6Convert(grid, basePars);
    dy = pchip_grad(grid, yGrid, x);
    dy(~isfinite(dy)) = 0;
    return
end
if iscell(L6pars{end}) && numel(L6pars{end}) >= 4 && strcmp(L6pars{end}{1}, 'postlocalbump')
    basePars = L6pars{end}{2};
    bumpPoints = L6pars{end}{3};
    bumpAmp = L6pars{end}{4};
    baseY = L6Convert(x, basePars);
    baseDy = L6Convert_grad(x, basePars);
    bumpY = local_bump_factor(x, bumpPoints, bumpAmp);
    bumpDy = local_bump_grad(x, bumpPoints, bumpAmp);
    dy = baseDy .* bumpY + baseY .* bumpDy;
    dy(~isfinite(dy)) = 0;
    return
end
if iscell(L6pars{end}) && numel(L6pars{end}) >= 4 && strcmp(L6pars{end}{1}, 'weightedblend')
    ratio = L6pars{end}{2};
    dy0 = L6Convert_grad(x, L6pars{end}{3});
    dy1 = L6Convert_grad(x, L6pars{end}{4});
    dy = (1-ratio) .* dy0 + ratio .* dy1;
    dy(~isfinite(dy)) = 0;
    return
end
h = 1e-5 * max(1, abs(x));
dy = (L6Convert(x + h, L6pars) - L6Convert(x - h, L6pars)) ./ (2*h);
dy(~isfinite(dy)) = 0;
end

function dy = pchip_grad(xGrid, yGrid, x)
[xGrid, sortIdx] = sort(xGrid(:).', 'ascend');
yGrid = yGrid(sortIdx);
pp = pchip(xGrid, yGrid);
breaks = pp.breaks;
coefs = pp.coefs;
nPieces = pp.pieces;
xFlat = x(:).';
dyFlat = zeros(size(xFlat));
for ii = 1:numel(xFlat)
    if xFlat(ii) <= breaks(1)
        piece = 1;
    elseif xFlat(ii) >= breaks(end)
        piece = nPieces;
    else
        piece = find(breaks <= xFlat(ii), 1, 'last');
        piece = min(piece, nPieces);
    end
    dx = xFlat(ii) - breaks(piece);
    c = coefs(piece, :);
    dyFlat(ii) = 3*c(1)*dx.^2 + 2*c(2)*dx + c(3);
end
dy = reshape(dyFlat, size(x));
end

function TweakFac = local_bump_factor(x,TweakPoints,TweakScale)
TweakFac = ones(size(x));
leftInd = x>=TweakPoints(1) & x<=TweakPoints(2);
rightInd = x>TweakPoints(2) & x<=TweakPoints(3);
TweakFac(leftInd) = 1 + TweakScale * 0.5 .* ...
    (1 - cos(pi*(x(leftInd) - TweakPoints(1))/(TweakPoints(2) - TweakPoints(1))));
TweakFac(rightInd) = 1 + TweakScale * 0.5 .* ...
    (1 + cos(pi*(x(rightInd) - TweakPoints(2))/(TweakPoints(3) - TweakPoints(2))));
end

function dTweak = local_bump_grad(x,TweakPoints,TweakScale)
dTweak = zeros(size(x));
leftInd = x>=TweakPoints(1) & x<=TweakPoints(2);
rightInd = x>TweakPoints(2) & x<=TweakPoints(3);
dTweak(leftInd) = TweakScale * 0.5 .* ...
    sin(pi*(x(leftInd) - TweakPoints(1))/(TweakPoints(2) - TweakPoints(1))) .* ...
    pi/(TweakPoints(2) - TweakPoints(1));
dTweak(rightInd) = -TweakScale * 0.5 .* ...
    sin(pi*(x(rightInd) - TweakPoints(2))/(TweakPoints(3) - TweakPoints(2))) .* ...
    pi/(TweakPoints(3) - TweakPoints(2));
end
