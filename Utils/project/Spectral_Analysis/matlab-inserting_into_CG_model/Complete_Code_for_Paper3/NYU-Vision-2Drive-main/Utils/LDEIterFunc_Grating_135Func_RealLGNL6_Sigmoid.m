%% This is the function to retreive books from the library
%% Here I use LocalResponse_3D to call Sigmoid functions corresponding to certain cells and LGN inputs

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

function LDEOutS = LDEIterFunc_Grating_135Func_RealLGNL6_Sigmoid(...
    celltype, L4EUse,L4IUse,...
    PixLGNCtgr,L6ELibInd)
    % Here I use LocalResponse_3D to find the corresponding Sigmoid
    % functions, the original version is to retrieve from library and do interpolations. 
    LDEOutLIBy = cell(5,1);
    LibyAll = zeros(length(L4EUse),5); % 5 LGN inputs (5 is the background LGN input)
    
    % 1-5 means 5 kinds of LGN inputs: oscillatory [90,0], [67.5,22.5], [45,45], [22.5,
    % 67.5]Hz, and Poisson 45Hz
    LDEOutLIBy{1} = ...
    LocalResponse_3D(celltype, L4EUse, L4IUse, 1, L6ELibInd);
    LibyAll(:,1) = LDEOutLIBy{1};
   
    LDEOutLIBy{2} = ...
    LocalResponse_3D(celltype, L4EUse, L4IUse, 2, L6ELibInd);
    LibyAll(:,2) = LDEOutLIBy{2};

    LDEOutLIBy{3} = ...
    LocalResponse_3D(celltype, L4EUse, L4IUse, 3, L6ELibInd);
    LibyAll(:,3) = LDEOutLIBy{3};

    LDEOutLIBy{4} = ...
    LocalResponse_3D(celltype, L4EUse, L4IUse, 4, L6ELibInd);
    LibyAll(:,4) = LDEOutLIBy{4};

    LDEOutLIBy{5} = ...
    LocalResponse_3D(celltype, L4EUse, L4IUse, 5, L6ELibInd);
    LibyAll(:,5) = LDEOutLIBy{5};
   
    % Then compose, to get the output firing rate. 
    LDEOutS = sum(PixLGNCtgr.*LibyAll,2);
    % 对不同的LGN做加权平均, 5个library 结果的加权平均值。

end