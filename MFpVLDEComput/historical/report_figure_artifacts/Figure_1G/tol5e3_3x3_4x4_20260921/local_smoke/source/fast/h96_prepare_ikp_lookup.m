function lookup = h96_prepare_ikp_lookup(IKp, step)
% Tabulate the one-dimensional IKp factor for fast uniform interpolation.
if nargin<2
    step = 0.001;
end
lookup.Min = 0;
lookup.Max = IKp.Highist+max(IKp.SmoothJoinHalfWidth)+1;
lookup.Step = step;
lookup.InvStep = 1/step;
lookup.Grid = (lookup.Min:step:lookup.Max).';
[~, lookup.Factor] = InhMulp(lookup.Grid, IKp);
lookup.LeftFactor = 1;
lookup.RightFactor = IKp.Slope;
end
