function output = h96_inh_mulp_fast(x, k)
% Active IKp7 multisigmoid, with constants and join polynomial precompiled.
firstZ = (x-k.Thrsld1)/(k.Thrsld2-k.Thrsld1)*(k.IntbH+k.IntbL)-k.IntbL;
firstProgress = (sigmoid(firstZ)-k.bLow)/(k.bHigh-k.bLow)*k.down1;
firstFactor = 1+firstProgress*(k.Slope-1);

secondZ = (x-k.Thrsld2)/(k.Highist-k.Thrsld2)*(k.IntaH+k.IntaL)-k.IntaL;
secondProgress = (sigmoid(secondZ)-k.aLow)/(k.aHigh-k.aLow)*k.down2+k.down1;
secondFactor = 1+secondProgress*(k.Slope-1);

halfWidths = k.SmoothJoinHalfWidth;
if isscalar(halfWidths)
    halfWidths = repmat(halfWidths, 1, 3);
end
factor = ones(size(x));
firstCore = x>=k.Thrsld1+halfWidths(1) & x<=k.Thrsld2-halfWidths(2);
secondCore = x>=k.Thrsld2+halfWidths(2) & x<=k.Highist-halfWidths(3);
tailCore = x>=k.Highist+halfWidths(3);
factor(firstCore) = firstFactor(firstCore);
factor(secondCore) = secondFactor(secondCore);
factor(tailCore) = k.Slope;

factor = blend(factor, firstFactor, x, k.Thrsld1, halfWidths(1));
if halfWidths(2)>0
    index = x>k.Thrsld2-halfWidths(2) & x<k.Thrsld2+halfWidths(2);
    if any(index, 'all')
        weight = compact_step((x(index)-(k.Thrsld2-halfWidths(2)))/(2*halfWidths(2)));
        factor(index) = (1-weight).*firstFactor(index)+weight.*secondFactor(index);
    end
end
tailFactor = k.Slope*ones(size(x));
index = x>k.Highist-halfWidths(3) & x<k.Highist+halfWidths(3);
if any(index, 'all')
    weight = compact_step((x(index)-(k.Highist-halfWidths(3)))/(2*halfWidths(3)));
    factor(index) = (1-weight).*secondFactor(index)+weight.*tailFactor(index);
end

index = x>=k.quinticLeft & x<=k.quinticRight;
if any(index, 'all')
    t = ((x(index)-k.quinticLeft)/(k.quinticRight-k.quinticLeft)).';
    powers = [ones(size(t)); t; t.^2; t.^3; t.^4; t.^5];
    factor(index) = sum(k.quinticCoeff.*powers, 1);
end
output = x.*factor;
end

function y = blend(y0, y1, x, center, halfWidth)
y = y0;
if halfWidth<=0
    return
end
index = x>center-halfWidth & x<center+halfWidth;
if any(index, 'all')
    weight = compact_step((x(index)-(center-halfWidth))/(2*halfWidth));
    y(index) = (1-weight).*y0(index)+weight.*y1(index);
end
end

function weight = compact_step(t)
weight = zeros(size(t));
weight(t>=1) = 1;
index = t>0 & t<1;
if any(index, 'all')
    a = exp(-1./t(index));
    b = exp(-1./(1-t(index)));
    weight(index) = a./(a+b);
end
end

function y = sigmoid(x)
y = exp(x)./(1+exp(x));
end
