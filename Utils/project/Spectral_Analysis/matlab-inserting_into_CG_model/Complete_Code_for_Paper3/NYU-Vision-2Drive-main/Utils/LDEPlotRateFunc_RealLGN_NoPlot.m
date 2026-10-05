function [LDEFrVec,ContourInfo] = LDEPlotRateFunc_RealLGN_NoPlot(...
    L4EPlot,L4IPlot,FrLDE,a1,a2,CtgrName,ODCtgr,CellCtgr,varargin)
if nargin > 8
    smth = varargin{1};
else
    smth = true;
end
if nargin >10
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
            CId1 = 1; CId2 = 2; WId = WIdS;
        case 'C'
            CId1 = 3; CId2 = 4; WId = WIdC;
        case 'I'
            CId1 = 5; CId2 = 5; WId = [0.5 0.5];
        otherwise
            error('No such cell category.')
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
            error('No such cell category.')
    end
end

if smth
    LDEFrcell = smoothdata(smoothdata(reshape(FrLDE(:,CId1)*WId(1)+FrLDE(:,CId2)*WId(2),a2,a1),2,'movmean',rowWinSize),'movmean',colWinSize);
else
    LDEFrcell = reshape(FrLDE(:,CId1)*WId(1)+FrLDE(:,CId2)*WId(2),a2,a1);
end

ContourInfo = [];
LDEFrVec = reshape(LDEFrcell,a1*a2,1);
end
