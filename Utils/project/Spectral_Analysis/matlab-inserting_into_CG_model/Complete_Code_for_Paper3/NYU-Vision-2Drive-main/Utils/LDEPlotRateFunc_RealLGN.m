%% Plot rate functions from LDE data
% Input: L4EPlot Grid for E Y*X
%        L4IPlot Grid for I Y*X
%        FrLDE   LDe computed firing rates (XY)*5: SOn COn SOff COff I
%        a1,a2   lengths of X Y range
%        CtgrName Name of OD Category
%        ODCtgr   OD Category: 1-3 
%        CellCtgr: S C I
%        varargin: smooth/not contour/not MFv/LIF [rowsmoothwin colsmoothwin]
% Output: LDEFrcell Fr function on grid L4EPlot,L4IPlot
%         ContourInfo: Linear slop/intersect of contour lines

function [LDEFrVec,ContourInfo] = LDEPlotRateFunc_RealLGN(...
    L4EPlot,L4IPlot,FrLDE,a1,a2,CtgrName,ODCtgr,CellCtgr,varargin)
if nargin > 8 % smooth or not
    smth = varargin{1};
else
    smth = true;
end

if nargin >9 % plot contour or not
    Contour = varargin{2};
else
    Contour = true;
end

if nargin >10 % realLGN (true, then 5 simulations) or fakeLGN (false, 3 simulations)
    realLGN = varargin{3};
else
    realLGN = true;
end

if nargin > 11
    rowWinSize = varargin{4}(1);
    colWinSize = varargin{4}(2);
else
    rowWinSize = 5;
    colWinSize = 5;
end

if realLGN
    Slgn = [4,6]; Clgn = [1,2]; 
    NSlgn = 4.5; NClgn = 1.5; 
    WIdS = (Slgn(2)-NSlgn)/(Slgn(2)-Slgn(1));
    WIdS = [WIdS,1-WIdS];
    WIdC = (Clgn(2)-NClgn)/(Clgn(2)-Clgn(1));
    WIdC = [WIdC,1-WIdC];
    switch CellCtgr
        case 'S'
            CId1 = 1; CId2 = 2;
            WId = WIdS;
        case 'C'
            CId1 = 3; CId2 = 4;
            WId = WIdC;
        case 'I'
            CId1 = 5; CId2 = 5;
            WId = [0.5 0.5];
        otherwise
            disp('No such cell category. Return')
            return
    end
else
    WId = [0.5 0.5];
    switch CellCtgr
        case 'S'
            CId1 = 1; CId2 = 1;
        case 'C'
            CId1 = 2; CId2 = 2;
        case 'I'
            CId1 = 3; CId2 = 3;
        otherwise
            disp('No such cell category. Return')
            return
    end
end

if smth
    LDEFrcell= smoothdata(...
                smoothdata(reshape(FrLDE(:,CId1)*WId(1)+FrLDE(:,CId2)*WId(2),a2,a1),...
                2,'movmean',rowWinSize),...
                  'movmean',colWinSize);
%  LDEFrcell= smoothdata(...
%                 smoothdata(reshape((FrLDE(:,CId1)+FrLDE(:,CId2))/2,a2,a1),2));
 else
    LDEFrcell =            reshape(FrLDE(:,CId1)*WId(1)+FrLDE(:,CId2)*WId(2),a2,a1);
end


% 2026/2/17 I commented the following three lines, I think we do not need
% them.
% s = mesh(L4EPlot,L4IPlot,LDEFrcell,'FaceAlpha','0.4');
% s.FaceColor = 'flat';
% s.EdgeColor = 'none';


hold on
if Contour
    [~,hh] = contour(L4EPlot,L4IPlot,LDEFrcell,5,'r',"ShowText",'on');
    view([0 90])
    xlabel('L4E'); ylabel('L4I');
    if iscell(CtgrName)
        title([CellCtgr ' ' CtgrName{ODCtgr}])
    elseif isnumeric(CtgrName)
        title(sprintf('%s LGN inpt %d, L6 inpt %d',CellCtgr,CtgrName,ODCtgr))
    end
    colorbar;
    axis tight
    
    % Export contour line slops
    ContourInfo = getContourLineCoordinates(hh);
%     ContourInfo = zeros(length(hh),2);
%     for ctInd = 1:length(hh)        
%         ctX = get(hh(ctInd),'XData');
%         ctY = get(hh(ctInd),'YData');
%         ContourInfo(ctInd,:) = polyfit(ctX,ctY,1);
%     end
else 
    ContourInfo = [];
end
LDEFrVec = reshape(LDEFrcell,a1*a2,1);
end