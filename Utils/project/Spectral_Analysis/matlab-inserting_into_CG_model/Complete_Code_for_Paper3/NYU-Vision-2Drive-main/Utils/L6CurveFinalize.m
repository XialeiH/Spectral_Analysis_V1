%% Finalized L6 curve. 
% xxLGN should be a vector where all entries in [0,1]

function [L6Out] = L6CurveFinalize(xxLGN)
L6Intesect = 0.65;
L6Shape = 0.18;
L6End = 1;
L6up = 60;
L6low = 6;

AngAbs = 0:7.5:180;
AngAbs = AngAbs - 90;
ScaledLGNL6_0deg = (cosd(abs(mod(AngAbs,180))*2)+1)/2;
L6_0deg = (rescaledSigmoid(ScaledLGNL6_0deg, L6Intesect, L6Shape) * L6End * (L6up-L6low) + L6low) /1e3;

AdjVec = [0.48,0.4,0.25,0.1,-0.4, -0.2, +0.05, -0.025];
AdjVecUse = [fliplr(AdjVec(2:end)),AdjVec];
L6_0New = L6_0deg;
Ind90 = floor(length(L6_0deg)/2)+1;
AdjInd = (Ind90-length(AdjVec)+1) : (Ind90+length(AdjVec)-1);
L6_0New(AdjInd) = L6_0New(AdjInd) + 0.06*0.3*AdjVecUse;

Y = fft(L6_0New(1:end-1));
N_dense = 1000; % For example, to make the output 5 times denser

% Prepare for higher resolution by zero-padding the frequency components
% The zeros are added symmetrically to maintain the signal's real-valued nature
Y_padded = [Y(1:floor(end/2)), zeros(1, N_dense-length(Y)), Y(floor(end/2)+1:end)];

% Perform IDFT on the modified frequency data
y_dense = ifft(Y_padded, 'symmetric');
y_dense = [y_dense,y_dense(1)];
y_dense = y_dense * max(L6_0New)/max(y_dense);
x_dense = linspace(-90,90,length(y_dense));

ScaledLGN = (cosd(abs(mod(x_dense,180))*2)+1)/2;
ScaledL6 = (y_dense * 1e3 - L6low)/(L6up - L6low);

ScaledLGNUse = sort(ScaledLGN(1:floor(N_dense/2)+1));
ScaledL6Use = sort(ScaledL6(1:floor(N_dense/2)+1));

%% finally, linearly strecthing the last part
L6LineInd = ScaledLGNUse>0.89;
Ind1 = find(L6LineInd);
Ind1 = Ind1(1)-1;
LineUp = max(ScaledL6Use);LineLow = ScaledL6Use(Ind1);
LineSlope = (LineUp - LineLow) / (ScaledLGNUse(end) - ScaledLGNUse(Ind1));
ScaledL6UseFix = ScaledL6Use;
ScaledL6UseFix(L6LineInd) = (ScaledLGNUse(L6LineInd) - ScaledLGNUse(Ind1)) * LineSlope + LineLow;
ScaledL6UseFix(ScaledL6UseFix<0) = 0;
L6Out = interp1(ScaledLGNUse, ScaledL6UseFix, xxLGN, 'Linear');
end

function y = rescaledSigmoid(x, xi, s)
    % Validate input
    if s < 0 || s > 1
        error('Shape parameter must be between 0 and 1');
    end
    if xi < 0 || xi > 1
        error('Intersection point must be between 0 and 1');
    end

    % Adjust k based on s
    k = 10 * (0.01 + s * 9);  % Scale factor for steepness

    % Calculate the value of L to ensure that the curve passes through (xi, xi)
    L = (1 / xi) - 1;

    % Apply the general sigmoid function
    y = L ./ (1 + exp(-k * (x - xi)));

    % Normalize to ensure the curve passes through (0,0) and (1,1)
    min_y = L / (1 + exp(k * xi));
    max_y = L / (1 + exp(-k * (1 - xi)));
    y = (y - min_y) / (max_y - min_y);
end

