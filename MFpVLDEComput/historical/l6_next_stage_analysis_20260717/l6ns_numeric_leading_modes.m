function modes = l6ns_numeric_leading_modes(phiFunction, fixedPoint, count, cfg)
% Matrix-free leading eigenmodes of D Phi at a nonlinear fixed point.

n = numel(fixedPoint);
base = phiFunction(fixedPoint);
step = 2e-6 * max(1, norm(fixedPoint) / sqrt(n));
operator = @(v) local_directional_action(phiFunction,fixedPoint,v,step,base);
opts = struct('tol',cfg.EigsTolerance,'maxit',cfg.EigsMaxIterations, ...
    'p',min(max(cfg.EigsSubspaceDimension,2*count+8),n), ...
    'issym',false,'isreal',true,'disp',0);
[vectors,values] = eigs(operator,n,count,'largestreal',opts);
lambda = diag(values);
[~,order] = sort(real(lambda),'descend');
modes = struct('Lambda',lambda(order),'Right',vectors(:,order));
end

function action=local_directional_action(phiFunction,fixedPoint,direction,step,base)
directionNorm=norm(direction);
if directionNorm==0
    action=zeros(size(base));
    return
end
h=step/directionNorm;
action=(phiFunction(fixedPoint+h*direction)-phiFunction(fixedPoint-h*direction))/(2*h);
end
