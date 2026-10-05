function y = l6ns_expmv_krylov(matrix, time, vector, krylovDimension, tolerance)
% Arnoldi approximation of exp(time*matrix)*vector.

if nargin < 4 || isempty(krylovDimension); krylovDimension = 60; end
if nargin < 5 || isempty(tolerance); tolerance = 1e-10; end

n = size(matrix, 1);
beta = norm(vector);
if beta == 0 || time == 0
    y = vector;
    return
end

m = min(krylovDimension, n);
v = zeros(n, m + 1, 'like', vector + 1i*0);
h = zeros(m + 1, m, 'like', vector + 1i*0);
v(:,1) = vector / beta;
used = m;

for j = 1:m
    candidate = matrix * v(:,j);
    for i = 1:j
        h(i,j) = v(:,i)' * candidate;
        candidate = candidate - h(i,j) * v(:,i);
    end
    h(j+1,j) = norm(candidate);
    if h(j+1,j) <= tolerance
        used = j;
        break
    end
    v(:,j+1) = candidate / h(j+1,j);
end

smallH = h(1:used, 1:used);
e1 = zeros(used, 1, 'like', smallH);
e1(1) = beta;
y = v(:,1:used) * (expm(time * smallH) * e1);
if isreal(matrix) && isreal(vector)
    y = real(y);
end
end
