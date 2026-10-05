function k = h96_prepare_ikp_fast(IKp)
% Compile the active smoothed multisigmoid IKp parameters to scalar constants.
assert(strcmp(IKp.Mode, 'multisigmoid'), ...
    'Fast IKp currently requires the active multisigmoid mode.');
k = IKp;
k.down2 = 1 - IKp.down1;
k.bHigh = sigmoid(IKp.IntbH);
k.bLow = sigmoid(-IKp.IntbL);
k.aHigh = sigmoid(IKp.IntaH);
k.aLow = sigmoid(-IKp.IntaL);

width = IKp.SmoothT2QuinticWidth;
k.quinticLeft = max(IKp.Thrsld1, IKp.Thrsld2 - width);
k.quinticRight = min(IKp.Highist, IKp.Thrsld2 + width);
dx = k.quinticRight - k.quinticLeft;
[y0, dy0, ddy0] = raw_derivatives(k.quinticLeft, k);
[y1, dy1, ddy1] = raw_derivatives(k.quinticRight, k);
c0 = y0;
c1 = dy0 * dx;
c2 = 0.5 * ddy0 * dx^2;
r1 = y1 - (c0 + c1 + c2);
r2 = dy1 * dx - (c1 + 2*c2);
r3 = ddy1 * dx^2 - 2*c2;
c5 = 0.5 * (r3 + 12*r1 - 6*r2);
c4 = 7*r2 - 15*r1 - r3;
c3 = 10*r1 - 4*r2 + 0.5*r3;
k.quinticCoeff = [c0; c1; c2; c3; c4; c5];
end

function [factor, slope, curvature] = raw_derivatives(x, k)
if x <= k.Thrsld2
    high = k.bHigh;
    low = k.bLow;
    dzdx = (k.IntbH+k.IntbL)/(k.Thrsld2-k.Thrsld1);
    z = (x-k.Thrsld1)*dzdx-k.IntbL;
    scale = k.down1/(high-low);
else
    high = k.aHigh;
    low = k.aLow;
    dzdx = (k.IntaH+k.IntaL)/(k.Highist-k.Thrsld2);
    z = (x-k.Thrsld2)*dzdx-k.IntaL;
    scale = k.down2/(high-low);
end
s = sigmoid(z);
progress = (s-low)*scale;
if x > k.Thrsld2
    progress = progress+k.down1;
end
slopeProgress = s*(1-s)*dzdx*scale;
curvatureProgress = s*(1-s)*(1-2*s)*dzdx^2*scale;
factor = 1+progress*(k.Slope-1);
slope = slopeProgress*(k.Slope-1);
curvature = curvatureProgress*(k.Slope-1);
end

function y = sigmoid(x)
y = exp(x)./(1+exp(x));
end
