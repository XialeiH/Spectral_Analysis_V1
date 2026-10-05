% L6fortyfive 1*N or N*1 vector!!
% MoveDeg: degree moved for the second curve, on a scale of 0 to 180 deg.

function L6zero = reverseL6Angle(L6fortyfive,MoveDeg)
N = length(L6fortyfive)-1;
%% L6fortyfive should be a 0:xx:180, so 1. get rid of 180, 2. get index
UnitDeg = 180/N;

MoveInd = floor(MoveDeg/UnitDeg);
% Construct matrix A
A = eye(N); % f(x) contribution
for i = 1:N
    prevIndex = i - MoveInd; % f(x - 45) index
    if prevIndex < 1
        prevIndex = prevIndex + N; % Wrap around for periodic boundary
    end
    A(i, prevIndex) = 1; % f(x - 45) contribution
end
%PurtA = 0.01*diag(eye(MoveInd)); % should be a row vec
PurtA = 0.05 * rand(1,MoveInd+1);

A(1:MoveInd+1,1:MoveInd+1) = A(1:MoveInd+1,1:MoveInd+1) + ...
    fliplr(eye(MoveInd+1)-diag(PurtA)) + diag(PurtA);

A = A/2; % since this is averaging two L6 domain

A_half = A(1:floor(N/2)+1, 1:floor(N/2)+1);


% Given G, solve for F
G = L6fortyfive; % Your g(x) vector here
%L6zeroPre = A\G(1:end-1)'; % Solve the linear system
L6zeroPre = lsqnonneg(A_half, G(1:floor(N/2)+1)');

L6zeroPre = [L6zeroPre(:);flipud(L6zeroPre(1:end-1))];
L6zero = reshape(L6zeroPre, size(L6fortyfive,1), size(L6fortyfive,2));
end