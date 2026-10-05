function d = HC_norm_diff(S1, C1, I1, S2, C2, I2, wC, aE, aI)
%HC_NORM_DIFF  HC-norm |f-g|_HC between two firing-rate configs.
%   S1,C1,I1 : state f   (one HC)
%   S2,C2,I2 : state g   (one HC)
%   wC       : weight of C in the effective E-population
%   aE,aI    : proportions of E- and I-cells in each pixel (aE+aI=1)

    % Effective E-populations
    E1 = (1 - wC) .* S1 + wC .* C1;
    E2 = (1 - wC) .* S2 + wC .* C2;

    % Differences per pixel (f - g)
    dE = E1 - E2;
    dI = I1 - I2;

    % Number of pixels in this HC
    Npix = numel(dE);

    % HC norm as defined in the paper
    d = sqrt( (1 / Npix) * sum( aE * dE(:).^2 + aI * dI(:).^2 ) );
end