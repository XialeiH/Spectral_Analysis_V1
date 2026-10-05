%% Iterations of LDE: Use precomputed function to determine output of a cell type
% Input: L4EmeshX, L4ImeshY  Domain of precomputed functions
%        LDEFrfunc_Subf      Functions of all input type for certain cell
%        population
%        L4EUse,L4IUse       All input L4E L4I
%        LDEUse              Just for output formality
%        varargin:
%        ORTPix,OBLPix,OPTPix,ORTOBLBd_Pix,OPTOBLBd_Pix   Pixel index of different ODs
%        or
%        PixLGNCtgr: n*3

% Zhuo-Cheng Xiao 04/06/2024

function LDEOutS = LDEIterFunc_Grating_135Func_RealLGNL6(...
    L4EmeshX, L4ImeshY,L6MeshZ,...
    LDEFrfunc_Subf,...
    L4EUse,L4IUse,...
    PixLGNCtgr,L6ELibInd,CtgrOrderReadout)
    % first look up the library
    LDEOutLIBy = cell(size(LDEFrfunc_Subf,1),1);
    LibyAll = zeros(length(L4EUse),size(LDEFrfunc_Subf,1));
    [XX,YY,ZZ] = meshgrid(unique(L4EmeshX), unique(L4ImeshY), L6MeshZ);
    for LGNInd = 1:size(LDEFrfunc_Subf,1)
        LDEOutLIBy{LGNInd} = ...
            interp3(XX,YY,ZZ,...
             squeeze(LDEFrfunc_Subf(LGNInd,:,:,:)),...
            L4EUse,   L4IUse,  L6ELibInd, 'linear');%,'makima'
        LibyAll(:,LGNInd) = LDEOutLIBy{LGNInd};
    end
    if numel(CtgrOrderReadout) ~= 4
        error('CtgrOrderReadout must contain 4 reordered foreground categories.');
    end

    % The library lookup first returns the canonical LGN1..4 firing rates.
    % Then reassign those four foreground channels according to the
    % preferred-angle / input-angle category mapping used in the
    % contr5/14 correction. Keep the background channel unchanged.
    LibyFG = LibyAll(:,1:4);
    LibyAll(:,1:4) = LibyFG(:, CtgrOrderReadout(:)');

    LDEOutS = sum(PixLGNCtgr.*LibyAll,2);




end
