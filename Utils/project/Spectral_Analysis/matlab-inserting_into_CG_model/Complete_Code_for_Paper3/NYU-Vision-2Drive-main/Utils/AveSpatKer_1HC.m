function C_SS_mean = AveSpatKer_1HC(C_PQ_Pixel_Us, N_HCin, N_HCout_in, NPixX, NPixY, varargin)
%% AveSpatKer_1HC
% Build an averaged spatial kernel, but ONLY for a single output HC.
% Output: C_SS_mean is (NPixX*NPixY) x (NPixX*NPixY).
%
% Input:
%   C_PQ_Pixel_Us : full pixel connectivity (Pre(col) * Post(row))
%   N_HCin        : # of HCs in the original network (e.g. 3 or 5)
%   N_HCout_in    : ignored (we force 1), kept for interface compatibility
%   NPixX, NPixY  : pixels per HC in x/y
%   varargin      : optional scaling constant of connectivity

%% Optional scale
if ~isempty(varargin)
    Scale = varargin{1};
else
    Scale = 1;  % default = no scaling
end

%% Basic checks
if size(C_PQ_Pixel_Us,1) ~= size(C_PQ_Pixel_Us,2)
    error('Input Matrix not square!')
end
if size(C_PQ_Pixel_Us,2) ~= N_HCin*NPixX*N_HCin*NPixY
    error('# of pixels doesnt match!')
end

if N_HCin < 3
    disp('Warning!: Reverse averaged kernel may have duplicate pixel infos');
end

if N_HCout_in ~= 1
    warning('AveSpatKer_1HC: forcing N_HCout = 1 (single HC output).');
end

%% 1) Build averaged 2HC-by-2HC kernel Ker_PQ_mean (same logic as original)

PrySynDist_all = zeros(2*NPixY, 2*NPixX, size(C_PQ_Pixel_Us,1));

for PixInd = 1:size(C_PQ_Pixel_Us,1)
    % spatial coord of this postsynaptic pixel in the original N_HCin x N_HCin HCs
    PX = ceil(PixInd/(N_HCin*NPixY)); 
    PY = mod(PixInd, N_HCin*NPixY);
    if PY == 0
        PY = N_HCin*NPixY;
    end
    
    % reshape whole map into [N_HCin*NPixY] x [N_HCin*NPixX]
    ConnVec = C_PQ_Pixel_Us(PixInd,:)';
    ConnMap = reshape(ConnVec, N_HCin*NPixY, N_HCin*NPixX);

    % center HC and bars
    Cen = floor((N_HCin-1)/2);
    HC_Center      = ConnMap(Cen*NPixY+1:(Cen+1)*NPixY, Cen*NPixX+1:(Cen+1)*NPixX);
    Bar_Center_hor = ConnMap(Cen*NPixY+1:(Cen+1)*NPixY, :);
    Bar_Center_ver = ConnMap(:, Cen*NPixX+1:(Cen+1)*NPixX);
    
    % extend by 1 HC on each side → (N_HCin+2) x (N_HCin+2) HC grid
    ConnMap_Ext = zeros((N_HCin+2)*NPixY, (N_HCin+2)*NPixX);
    ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY, ...
                NPixX+1:(N_HCin+1)*NPixX) = ConnMap;
    % 4 corners = center HC
    ConnMap_Ext(1:NPixY, 1:NPixX) = HC_Center;
    ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY, 1:NPixX) = HC_Center;
    ConnMap_Ext(1:NPixY, (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = HC_Center;
    ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY, ...
                (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = HC_Center;
    % 4 sides = center bars
    ConnMap_Ext(1:NPixY, NPixX+1:(N_HCin+1)*NPixX) = Bar_Center_hor;
    ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY, ...
                NPixX+1:(N_HCin+1)*NPixX) = Bar_Center_hor;
    ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY, 1:NPixX) = Bar_Center_ver;
    ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY, ...
                (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = Bar_Center_ver;
    
    % Extract the 2HC×2HC patch around this pixel
    PX_New = PX + NPixX;
    PY_New = PY + NPixY;
    PrySynDist_all(:,:,PixInd) = ConnMap_Ext(PY_New-NPixY:PY_New+NPixY-1, ...
                                             PX_New-NPixX:PX_New+NPixX-1);
end

Ker_PQ_mean = mean(PrySynDist_all, 3);

%% 2) Symmetrize Ker_PQ_mean (FIXED: allocate 4 slices for diagonal/antidiagonal)

% flip up-down
Ker_PQ_meanSym = zeros(2*NPixY, 2*NPixX, 2);
Ker_PQ_meanSym(2:end,2:end,1) = Ker_PQ_mean(2:end,2:end);
Ker_PQ_meanSym(2:end,2:end,2) = Ker_PQ_mean(end:-1:2,2:end);
Ker_PQ_mean = mean(Ker_PQ_meanSym,3);

% flip left-right
Ker_PQ_meanSym = zeros(2*NPixY,2*NPixX,2);
Ker_PQ_meanSym(2:end,2:end,1) = Ker_PQ_mean(2:end,end:-1:2);
Ker_PQ_meanSym(2:end,2:end,2) = Ker_PQ_mean(2:end,2:end);
Ker_PQ_mean = mean(Ker_PQ_meanSym,3);

% diagonal + antidiagonal symmetrization
Ker_PQ_meanDia = zeros(2*NPixY,2*NPixX,4);   % <-- must be 4, not 2
Ker_PQ_meanDia(2:end,2:end,1) = Ker_PQ_mean(2:end,2:end);
Ker_PQ_meanDia(2:end,2:end,2) = Ker_PQ_mean(2:end,2:end)';                  % transpose
Ker_PQ_meanDia(2:end,2:end,3) = rot90(Ker_PQ_mean(2:end,2:end),2);          % 180° rotation
Ker_PQ_meanDia(2:end,2:end,4) = Ker_PQ_mean(2:end,2:end);
Ker_PQ_mean = mean(Ker_PQ_meanDia,3);

% apply global scale
Ker_PQ_mean = Scale * Ker_PQ_mean;

%% 3) Put kernel back to matrix, ONLY for 1 HC with periodic BC
%    IMPORTANT FIX: use *summing* under periodic wrapping, not overwriting.

N_HCout  = 1;                     % force 1-HC output conceptually
Npix_out = NPixX * NPixY;         % #pixels in the central HC

C_SS_mean = zeros(Npix_out, Npix_out);   % rows = post, cols = pre

for PixInd = 1:Npix_out
    % Pixel coordinates (PX, PY) inside the single HC:
    % X runs 1..NPixX, Y runs 1..NPixY
    PX = ceil(PixInd / NPixY);     % column index within 1 HC
    PY = mod(PixInd, NPixY);       % row index within 1 HC
    if PY == 0
        PY = NPixY;
    end

    % Periodic embedding of the 2HC×2HC kernel Ker_PQ_mean onto 1HC torus
    ConnMap_Rev = zeros(NPixY, NPixX);

    % Loop over the 2NPixY x 2NPixX kernel and wrap onto the 1HC lattice
    for dy = -NPixY:(NPixY-1)
        ky = dy + NPixY + 1;    % 1..2NPixY
        for dx = -NPixX:(NPixX-1)
            kx = dx + NPixX + 1;   % 1..2NPixX

            % periodic y/x indices in the 1-HC domain
            yy = mod(PY + dy - 1, NPixY) + 1;
            xx = mod(PX + dx - 1, NPixX) + 1;

            % ADD contribution (not overwrite!)
            ConnMap_Rev(yy, xx) = ConnMap_Rev(yy, xx) + Ker_PQ_mean(ky, kx);
        end
    end

    % Flatten to one row of the final matrix (pre-synaptic index along columns)
    C_SS_mean(PixInd, :) = reshape(ConnMap_Rev, 1, Npix_out);
end

end










%% commented below is the old version with some problems
% function C_SS_mean = AveSpatKer_1HC(C_PQ_Pixel_Us, N_HCin, N_HCout_in, NPixX, NPixY, varargin)
% %% AveSpatKer_1HC
% % Build an averaged spatial kernel, but ONLY for a single output HC.
% % Output: C_SS_mean is (NPixX*NPixY) x (NPixX*NPixY).
% %
% % Input:
% %   C_PQ_Pixel_Us : full pixel connectivity (Pre(col) * Post(row))
% %   N_HCin        : # of HCs in the original network (e.g. 3)
% %   N_HCout_in    : ignored, but kept for interface compatibility
% %   NPixX, NPixY  : pixels per HC in x/y
% %   varargin      : scaling constant of connectivity
% 
% if ~isempty(varargin)
%     Scale = varargin{1};
% else
%     Scale = 1;  % default = no scaling
% end
% 
% %% Basic checks
% if size(C_PQ_Pixel_Us,1) ~= size(C_PQ_Pixel_Us,2)
%     error('Input Matrix not square!')
% end
% if size(C_PQ_Pixel_Us,2) ~= N_HCin*NPixX*N_HCin*NPixY
%     error('# of pixels doesnt match!')
% end
% 
% if N_HCin < 3
%     disp('Warning!: Reverse averaged kernel may have duplicate pixel infos');
% end
% 
% if N_HCout_in ~= 1
%     warning('AveSpatKer_1HC: forcing N_HCout = 1 (single HC output).');
% end
% 
% %% 1) Build averaged 2HC-by-2HC kernel Ker_PQ_mean (same as original)
% 
% PrySynDist_all = zeros(2*NPixY, 2*NPixX, size(C_PQ_Pixel_Us,1));
% 
% for PixInd = 1:size(C_PQ_Pixel_Us,1)
%     % spatial coord of this postsynaptic pixel in the original N_HCin x N_HCin HCs
%     PX = ceil(PixInd/(N_HCin*NPixY)); 
%     PY = mod(PixInd, N_HCin*NPixY);
%     if PY == 0
%         PY = N_HCin*NPixY;
%     end
% 
%     % reshape whole map
%     ConnVec = C_PQ_Pixel_Us(PixInd,:)';
%     ConnMap = reshape(ConnVec, N_HCin*NPixY, N_HCin*NPixX);
% 
%     % center HC and bars
%     Cen = floor((N_HCin-1)/2);
%     HC_Center      = ConnMap(Cen*NPixY+1:(Cen+1)*NPixY, Cen*NPixX+1:(Cen+1)*NPixX);
%     Bar_Center_hor = ConnMap(Cen*NPixY+1:(Cen+1)*NPixY, :);
%     Bar_Center_ver = ConnMap(:, Cen*NPixX+1:(Cen+1)*NPixX);
% 
%     % extend by 1 HC on each side → (N_HCin+2) x (N_HCin+2) HC grid
%     ConnMap_Ext = zeros((N_HCin+2)*NPixY, (N_HCin+2)*NPixX);
%     ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY, ...
%                 NPixX+1:(N_HCin+1)*NPixX) = ConnMap;
%     % 4 corners = center HC
%     ConnMap_Ext(1:NPixY, 1:NPixX) = HC_Center;
%     ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY, 1:NPixX) = HC_Center;
%     ConnMap_Ext(1:NPixY, (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = HC_Center;
%     ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY, ...
%                 (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = HC_Center;
%     % 4 sides = center bars
%     ConnMap_Ext(1:NPixY, NPixX+1:(N_HCin+1)*NPixX) = Bar_Center_hor;
%     ConnMap_Ext((N_HCin+1)*NPixY+1:(N_HCin+2)*NPixY, ...
%                 NPixX+1:(N_HCin+1)*NPixX) = Bar_Center_hor;
%     ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY, 1:NPixX) = Bar_Center_ver;
%     ConnMap_Ext(NPixY+1:(N_HCin+1)*NPixY, ...
%                 (N_HCin+1)*NPixX+1:(N_HCin+2)*NPixX) = Bar_Center_ver;
% 
%     % Extract the 2HC×2HC patch around this pixel
%     PX_New = PX + NPixX;
%     PY_New = PY + NPixY;
%     PrySynDist_all(:,:,PixInd) = ConnMap_Ext(PY_New-NPixY:PY_New+NPixY-1, ...
%                                              PX_New-NPixX:PX_New+NPixX-1);
% end
% 
% Ker_PQ_mean = mean(PrySynDist_all, 3);
% 
% %% 2) Symmetrize Ker_PQ_mean (same as original)
% 
% Ker_PQ_meanSym = zeros(2*NPixY, 2*NPixX, 2);
% Ker_PQ_meanSym(2:end,2:end,1) = Ker_PQ_mean(2:end,2:end);
% Ker_PQ_meanSym(2:end,2:end,2) = Ker_PQ_mean(end:-1:2,2:end);
% Ker_PQ_mean = mean(Ker_PQ_meanSym,3);
% 
% Ker_PQ_meanSym = zeros(2*NPixY,2*NPixX,2);
% Ker_PQ_meanSym(2:end,2:end,1) = Ker_PQ_mean(2:end,end:-1:2);
% Ker_PQ_meanSym(2:end,2:end,2) = Ker_PQ_mean(2:end,2:end);
% Ker_PQ_mean = mean(Ker_PQ_meanSym,3);
% 
% Ker_PQ_meanDia = zeros(2*NPixY,2*NPixX,2);
% Ker_PQ_meanDia(2:end,2:end,1) = Ker_PQ_mean(2:end,2:end);
% Ker_PQ_meanDia(2:end,2:end,2) = Ker_PQ_mean(2:end,2:end)';
% Ker_PQ_mean = mean(Ker_PQ_meanDia,3);
% 
% Ker_PQ_meanDia(2:end,2:end,3) = rot90(Ker_PQ_mean(2:end,2:end),2);
% Ker_PQ_meanDia(2:end,2:end,4) = Ker_PQ_mean(2:end,2:end);
% Ker_PQ_mean = mean(Ker_PQ_meanDia,3);
% 
% Ker_PQ_mean = Scale * Ker_PQ_mean;
% 
% %% 3) Put kernel back to matrix, but ONLY for 1 HC with periodic BC
% 
% N_HCout  = 1;                     % force 1-HC output conceptually
% Npix_out = NPixX * NPixY;         % #pixels in the central HC
% 
% C_SS_mean = zeros(Npix_out, Npix_out);   % Pre(col) × Post(row) or row×col depending on your convention
% 
% for PixInd = 1:Npix_out
%     % Pixel coordinates (PX, PY) inside the single HC:
%     % X runs 1..NPixX, Y runs 1..NPixY
%     PX = ceil(PixInd / NPixY);     % column index within 1 HC
%     PY = mod(PixInd, NPixY);       % row index within 1 HC
%     if PY == 0
%         PY = NPixY;
%     end
% 
%     % Periodic embedding of the 2HC×2HC kernel Ker_PQ_mean onto 1HC torus
%     ConnMap_Rev = zeros(NPixY, NPixX);
% 
%     % Indices for the 2HC×2HC window wrapped back to 1HC (periodic BC)
%     YRange = mod(PY-NPixY:PY+NPixY-1, NPixY);
%     YRange(YRange == 0) = NPixY;
%     XRange = mod(PX-NPixX:PX+NPixX-1, NPixX);
%     XRange(XRange == 0) = NPixX;
% 
%     % Place the 20×20 kernel onto the 10×10 HC with wrapping
%     ConnMap_Rev(YRange, XRange) = Ker_PQ_mean;
% 
%     % Flatten to one row of the final matrix
%     C_SS_mean(PixInd, :) = reshape(ConnMap_Rev, 1, Npix_out);
% end
% 
% end
