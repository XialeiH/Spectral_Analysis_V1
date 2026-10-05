function LDE_1HC = extract_left_optimal_HC(LDE_equiv, N_HCin, NPixX, NPixY)
% LDE_equiv: struct with fields S, C, I, each length (N_HCin^2 * NPixX * NPixY)
% N_HCin   : number of HCs in original data (e.g. 4)
% NPixX/Y  : pixels per HC in x/y
%
% LDE_1HC  : struct with fields S, C, I, each length (NPixX * NPixY),
%            containing only the central HC.

    Nx = N_HCin * NPixX;
    Ny = N_HCin * NPixY;

    % "Central" HC index in each direction, same convention as AveSpatKer:
    % Cen = floor((N_HCin-1)/2).
    Cen = floor((N_HCin-1)/2);  % for N_HCin=4 -> Cen=1 -> HC #2

    x_idx = (Cen+1)*NPixX + (1:NPixX);   % columns of central HC
    y_idx = Cen*NPixY + (1:NPixY);   % rows of central HC

    fields = {'S','C','I'};
    for k = 1:numel(fields)
        Q = fields{k};
        full_vec = LDE_equiv.(Q);              % length N_HCin^2 * NPixX * NPixY
        full_map = reshape(full_vec, Ny, Nx); % Ny × Nx (e.g. 40×40)

        block = full_map(y_idx, x_idx);       % central HC (10×10)
        LDE_1HC.(Q) = block(:);               % flatten to 100×1
    end
end
