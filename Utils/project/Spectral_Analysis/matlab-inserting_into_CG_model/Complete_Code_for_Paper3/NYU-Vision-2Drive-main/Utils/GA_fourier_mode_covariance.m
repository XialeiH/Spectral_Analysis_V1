function SigmaHat = GA_fourier_mode_covariance(Ahat, Qhat, mode)
% Mode-wise covariance solve for approximately translation-invariant sheets.
% Ahat/Qhat can be cell arrays over sheet modes or small 3D arrays.
if nargin < 3 || isempty(mode)
    mode = 'discrete';
end

if iscell(Ahat)
    SigmaHat = cell(size(Ahat));
    for k = 1:numel(Ahat)
        SigmaHat{k} = GA_linearized_covariance(Ahat{k}, Qhat{k}, mode);
    end
else
    nMode = size(Ahat, 3);
    SigmaHat = zeros(size(Ahat));
    for k = 1:nMode
        SigmaHat(:,:,k) = GA_linearized_covariance(Ahat(:,:,k), Qhat(:,:,k), mode);
    end
end
end
