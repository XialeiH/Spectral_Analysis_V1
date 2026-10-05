function validation = l6ns_true_parameter_validation(cfg,data,context)
% Sparse moved-fixed-point validation of the two FPP pathway predictions.

outputDir = fullfile(cfg.OutputRoot,'true_parameter_validation');
if ~exist(outputDir,'dir'); mkdir(outputDir); end
pathway = l6ns_two_pathway_data(data);
fixed = context.FixedPoint(:);
rawContext = context;
rawContext.FixedPointCorrection = [];
context.FixedPointCorrection = fixed-l6ns_phi(fixed,0,rawContext,[1 1],1);

pathwayNames = {'L6','I'};
gainValues = {[0.95 1.05 1.15],[0.95 1.05 0.90]};
rows = {};
for pathwayIndex = 1:2
    pathwayName = pathwayNames{pathwayIndex};
    oneGainGrid = gainValues{pathwayIndex};
    for gainIndex = 1:numel(oneGainGrid)
        gain = oneGainGrid(gainIndex);
        if strcmp(pathwayName,'L6')
            phiTrue = @(x)l6ns_phi(x,1-gain,context,[1 1],1,'true');
            jFpp = l6ns_pathway_jacobian(pathway,gain,1);
        else
            phiTrue = @(x)l6ns_phi(x,0,context,[1 1],gain,'true');
            jFpp = l6ns_pathway_jacobian(pathway,1,gain);
        end
        [moved,residual,iterations] = local_fixed_point(phiTrue,fixed,context.RelaxationP);
        [lambdaTrue,rightTrue,eigenResidual] = local_matrix_free_leading(phiTrue,moved);
        modesFpp = l6ns_eigenpairs(jFpp,6,'largestreal',cfg);
        [lambdaFpp,indexFpp] = max(real(modesFpp.Lambda));
        rightFpp = modesFpp.Right(:,indexFpp);
        overlap = abs(rightTrue'*rightFpp)/max(norm(rightTrue)*norm(rightFpp),eps);
        rows(end+1,:) = {pathwayName,gain,norm(moved-fixed)/max(norm(fixed),eps), ...
            residual,iterations,real(lambdaTrue),imag(lambdaTrue),lambdaFpp, ...
            real(lambdaTrue)-lambdaFpp,overlap,eigenResidual}; %#ok<AGROW>
    end
end
summary = cell2table(rows,'VariableNames', ...
    {'pathway','gain','relativeFixedPointShift','fixedPointResidual','iterations', ...
    'trueLambdaReal','trueLambdaImag','fppLambdaReal','lambdaDifference', ...
    'rightModeOverlap','matrixFreeEigenResidual'});
writetable(summary,fullfile(outputDir,'true_vs_fpp_validation.tsv'), ...
    'FileType','text','Delimiter','\t');

fig = figure('Visible','off','Color','w','Position',[100 100 1000 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for pathwayIndex = 1:2
    nexttile;
    one = summary(strcmp(summary.pathway,pathwayNames{pathwayIndex}),:);
    one = sortrows(one,'gain');
    plot(one.gain,one.trueLambdaReal,'-o','LineWidth',1.5); hold on;
    plot(one.gain,one.fppLambdaReal,'--o','LineWidth',1.5); grid on;
    xlabel('pathway gain'); ylabel('max Re \lambda');
    title(sprintf('%s true versus FPP',pathwayNames{pathwayIndex}));
    legend({'true moved fixed point','FPP'},'Location','best');
end
l6ns_save_figure(fig,outputDir,'figure6_true_vs_fpp_validation'); close(fig);
validation = struct('Summary',summary);
save(fullfile(outputDir,'true_parameter_validation_result.mat'),'validation','-v7.3');
end

function [state,residual,iteration] = local_fixed_point(phi,initial,p)
state = initial;
for iteration = 1:250
    next = (1-p)*state+p*phi(state);
    residual = norm(phi(next)-next)/max(norm(next),eps);
    state = next;
    if residual < 1e-10
        break
    end
end
end

function [lambda,right,residual] = local_matrix_free_leading(phi,state)
n = numel(state);
phiState = phi(state);
stateScale = max(1,norm(state)/sqrt(n));
operator = @(v)local_jvp(phi,state,phiState,v,stateScale);
opts = struct('tol',1e-6,'maxit',350,'p',40,'disp',0);
[vectors,values] = eigs(operator,n,4,'largestreal',opts);
lambdaPool = diag(values);
[~,index] = max(real(lambdaPool));
lambda = lambdaPool(index);
right = vectors(:,index); right=right/norm(right);
jRight = operator(right);
residual = norm(jRight-lambda*right)/max(norm(jRight),eps);
if residual > 5e-4
    warning('L6NS:MatrixFreeResidual','Matrix-free eigen residual is %.3e.',residual);
end
end

function output = local_jvp(phi,state,phiState,direction,stateScale)
if ~isreal(direction)
    output = local_jvp(phi,state,phiState,real(direction),stateScale) + ...
        1i*local_jvp(phi,state,phiState,imag(direction),stateScale);
    return
end
directionNorm = norm(direction);
if directionNorm == 0
    output = zeros(size(direction));
    return
end
h = 2e-6*stateScale/directionNorm;
output = (phi(state+h*direction)-phiState)/h;
end
