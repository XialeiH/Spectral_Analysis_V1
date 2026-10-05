


function C_SS_mean = AveSpatKer_1HC_Flip(C_PQ_Pixel_Us, N_HCin, N_HCout_in, NPixX, NPixY, all_connectivity_scale, single_pixel_connection_scale)



Scale = all_connectivity_scale;  

single_pixel_scale = single_pixel_connection_scale;


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

% Total pixels and basic sizes
Npix_total = size(C_PQ_Pixel_Us, 1);
Nrow_full  = N_HCin * NPixY;
Ncol_full  = N_HCin * NPixX;

% Index range of the central HC in the full [Nrow_full x Ncol_full] map
Cen         = floor((N_HCin-1)/2);
x_cen_start = Cen * NPixX + 1;
x_cen_end   = (Cen+1) * NPixX;
y_cen_start = Cen * NPixY + 1;
y_cen_end   = (Cen+1) * NPixY;

% We now only average over pixels in the *central* HC
Npix_central = NPixX * NPixY;
PrySynDist_all = zeros(2*NPixY, 2*NPixX, Npix_central);

idx_cen = 0;

for PixInd = 1:Npix_total
    % spatial coord of this postsynaptic pixel in the original N_HCin x N_HCin HCs
    PX = ceil(PixInd / Nrow_full);  % column index (x)
    PY = mod(PixInd, Nrow_full);    % row index (y)
    if PY == 0
        PY = Nrow_full;
    end

    % Only keep pixels that lie in the central HC
    if PX < x_cen_start || PX > x_cen_end || PY < y_cen_start || PY > y_cen_end
        continue
    end

    idx_cen = idx_cen + 1;

    % reshape whole map into [Nrow_full] x [Ncol_full]
    ConnVec = C_PQ_Pixel_Us(PixInd, :)';
    ConnMap = reshape(ConnVec, Nrow_full, Ncol_full);

    % Build 2HC×2HC patch around (PX,PY) *without* any extended tiling.
    % Outside the [1..Nrow_full]×[1..Ncol_full] region → zero (already default).
    patch = zeros(2*NPixY, 2*NPixX);

    % global window in the full map
    x_win = (PX - NPixX + 1):(PX + NPixX);
    y_win = (PY - NPixY + 1):(PY + NPixY);

    % intersect with valid indices in the full map
    x_valid = max(x_win(1), 1) : min(x_win(end), Ncol_full);
    y_valid = max(y_win(1), 1) : min(y_win(end), Nrow_full);

    % corresponding indices in the local 2HC×2HC patch
    x_patch = (x_valid - x_win(1)) + 1;
    y_patch = (y_valid - y_win(1)) + 1;

    % fill only the overlapping region; the rest stays zero
    patch(y_patch, x_patch) = ConnMap(y_valid, x_valid);

    PrySynDist_all(:, :, idx_cen) = patch;
end

Ker_PQ_mean = mean(PrySynDist_all, 3);
Ker_PQ_mean = Scale * Ker_PQ_mean;  % apply optional scaling if desired


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


Ker_PQ_mean



%% Build equivalent 1-HC connectivity by folding 9 HCs into the central HC

PixNum1HC = NPixX * NPixY;
C_SS_mean     = zeros(PixNum1HC, PixNum1HC);   % pre (row) × post (col) or vice versa as you like

% center of the 2HC×2HC kernel in pixel coordinates
cx = NPixX;
cy = NPixY;

% --- loop over postsynaptic pixels in the central HC --- %
for x_post = 1:NPixX
    for y_post = 1:NPixY

        % 1) position of this postsyn pixel in the 3×3 HC layout
        % central HC is at (HCx,HCy) = (2,2)
        x0_glob = NPixX + x_post;   % global x index (1..3*NPixX)
        y0_glob = NPixY + y_post;   % global y index (1..3*NPixY)

        % one-HC map to accumulate folded connections for this post pixel
        Conn1HC = zeros(NPixY, NPixX);

        % 2) loop over presynaptic pixels in the 2HC×2HC spatial kernel
        for dx = -(NPixX-1):NPixX
            for dy = -(NPixY-1):NPixY

                % location in the kernel
                kx = cx + dx;
                ky = cy + dy;

                % skip if outside kernel array (should not happen)
                if kx < 1 || kx > 2*NPixX || ky < 1 || ky > 2*NPixY
                    continue
                end

                w = Ker_PQ_mean(ky, kx);  % connection weight
                if w == 0
                    continue
                end

                % global presynaptic pixel index in the 3×3 HCs
                x_pre_glob = x0_glob + dx;
                y_pre_glob = y0_glob + dy;

                % which HC does this presynaptic pixel lie in? (1..3)
                HCx = ceil(x_pre_glob / NPixX);
                HCy = ceil(y_pre_glob / NPixY);

                % local coordinates inside that HC (1..NPixX/NPixY)
                x_in = x_pre_glob - (HCx-1)*NPixX;
                y_in = y_pre_glob - (HCy-1)*NPixY;

                % 3 & 4) move / flip presynaptic pixel into the central HC
                if HCx == 2 && HCy == 2
                    % central HC: keep as is
                    x_fold = x_in;
                    y_fold = y_in;

                elseif HCx == 1 && HCy == 2
                    % left edge: flip horizontally
                    x_fold = NPixX + 1 - x_in;
                    y_fold = y_in;

                elseif HCx == 3 && HCy == 2
                    % right edge: flip horizontally
                    x_fold = NPixX + 1 - x_in;
                    y_fold = y_in;

                elseif HCx == 2 && HCy == 1
                    % top edge: flip vertically
                    x_fold = x_in;
                    y_fold = NPixY + 1 - y_in;

                elseif HCx == 2 && HCy == 3
                    % bottom edge: flip vertically
                    x_fold = x_in;
                    y_fold = NPixY + 1 - y_in;

                elseif HCx == 1 && HCy == 1
                    % top-left corner: flip in both x and y (point reflection)
                    x_fold = NPixX + 1 - x_in;
                    y_fold = NPixY + 1 - y_in;

                elseif HCx == 3 && HCy == 1
                    % top-right corner
                    x_fold = NPixX + 1 - x_in;
                    y_fold = NPixY + 1 - y_in;

                elseif HCx == 1 && HCy == 3
                    % bottom-left corner
                    x_fold = NPixX + 1 - x_in;
                    y_fold = NPixY + 1 - y_in;

                elseif HCx == 3 && HCy == 3
                    % bottom-right corner
                    x_fold = NPixX + 1 - x_in;
                    y_fold = NPixY + 1 - y_in;

                else
                    % should not happen for a 2HC×2HC kernel on a 3×3 HC grid
                    continue
                end

                % 5) accumulate connection weight into the common central HC map
                Conn1HC(y_fold, x_fold) = Conn1HC(y_fold, x_fold) + w;

            end
        end

        % 6) reshape this 1HC map into one column of the final C
        post_ind = (x_post-1)*NPixY + y_post;      % column index for this post pixel
        C_SS_mean(:, post_ind) = Conn1HC(:);           % vectorized in (y,x) order

    end
end