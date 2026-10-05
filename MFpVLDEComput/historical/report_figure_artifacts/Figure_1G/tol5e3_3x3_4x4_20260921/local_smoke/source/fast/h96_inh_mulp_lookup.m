function output = h96_inh_mulp_lookup(x, lookup)
% Apply a pretabulated IKp factor with uniform linear interpolation.
factor = lookup.LeftFactor*ones(size(x));
factor(x>=lookup.Max) = lookup.RightFactor;
index = x>lookup.Min & x<lookup.Max;
if any(index, 'all')
    u = (x(index)-lookup.Min).*lookup.InvStep;
    offset = floor(u);
    fraction = u-offset;
    left = offset+1;
    factor(index) = lookup.Factor(left).*(1-fraction)+lookup.Factor(left+1).*fraction;
end
output = x.*factor;
end
