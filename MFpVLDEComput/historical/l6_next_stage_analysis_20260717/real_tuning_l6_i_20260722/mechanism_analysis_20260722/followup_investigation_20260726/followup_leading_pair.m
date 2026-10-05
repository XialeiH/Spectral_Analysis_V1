function [lambda,right,left,residual] = followup_leading_pair(J)
% Return a normalized leading left/right eigenpair of a sparse map Jacobian.

options = struct('tol',1e-10,'maxit',2600,'p',120,'isreal',true,'disp',0);
[rightVectors,rightValues] = eigs(J,16,'largestreal',options);
rightValues = diag(rightValues);
[~,index] = max(real(rightValues));
lambda = rightValues(index);
right = rightVectors(:,index);
[leftVectors,leftValues] = eigs(J',20,'largestreal',options);
leftValues = diag(leftValues);
[~,leftIndex] = min(abs(leftValues-conj(lambda)));
left = leftVectors(:,leftIndex);
normalizer = left'*right;
left = left/conj(normalizer);
residual = norm(J*right-lambda*right)/max(norm(J*right),eps);
end
