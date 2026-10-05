function Sigma = GA_linearized_covariance(A, Q, mode)
% Linearized response covariance for proposal section 5.6.
% mode='discrete': Sigma = A*Sigma*A' + Q.
% mode='continuous': (A-I)*Sigma + Sigma*(A-I)' + Q = 0.
if nargin < 3 || isempty(mode)
    mode = 'discrete';
end

switch lower(mode)
    case 'discrete'
        if exist('dlyap', 'file') == 2
            Sigma = dlyap(full(A), full(Q));
        else
            n = size(A, 1);
            Sigma = reshape((eye(n*n) - kron(full(A), full(A))) \ Q(:), n, n);
        end
    case 'continuous'
        L = full(A) - eye(size(A, 1));
        if exist('lyap', 'file') == 2
            Sigma = lyap(L, full(Q));
        else
            n = size(A, 1);
            K = kron(eye(n), L) + kron(L, eye(n));
            Sigma = reshape(-(K \ Q(:)), n, n);
        end
    otherwise
        error('Unknown covariance mode %s.', mode);
end

Sigma = 0.5 * (Sigma + Sigma');
end
