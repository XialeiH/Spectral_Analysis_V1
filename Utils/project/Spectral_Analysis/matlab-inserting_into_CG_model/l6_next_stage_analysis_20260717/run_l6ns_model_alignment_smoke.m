% Verify that the nonlinear Phi evaluator matches the saved endpoint Jacobians.

cfg=l6ns_config();
data=l6ns_load_endpoints(cfg);
canonicalScript=getenv('L6NS_CANONICAL_MODEL_SCRIPT');
if isempty(canonicalScript); error('L6NS_CANONICAL_MODEL_SCRIPT is required.'); end
if isempty(getenv('L6NS_CANONICAL_CODE_DIR')); error('L6NS_CANONICAL_CODE_DIR is required.'); end
setenv('GEOM_ANGLE_COND',sprintf('%.12g',cfg.Angle));
setenv('GEOM_CONTR_COND',sprintf('%.12g',cfg.Contrast));
setenv('L6_EQU_BLEND_WEIGHT','0');
setenv('L6NS_SETUP_ONLY','1');
run(canonicalScript);
run(fullfile(cfg.CodeRoot,'l6ns_context_from_canonical_workspace.m'));

fixed=ModelContext.FixedPoint;
weights=[0 1];
jacobians={data.J0,data.A};
rows=cell(2,1);
rng(cfg.RandomSeed);
for endpoint=1:2
    w=weights(endpoint);
    phi=@(x)l6ns_phi(x,w,ModelContext);
    fixedResidual=norm(phi(fixed)-fixed)/max(norm(fixed),eps);
    errors=nan(4,1);
    step=2e-6*max(1,norm(fixed)/sqrt(numel(fixed)));
    for directionIndex=1:numel(errors)
        direction=randn(size(fixed)); direction=direction/norm(direction);
        numeric=(phi(fixed+step*direction)-phi(fixed-step*direction))/(2*step);
        analytic=jacobians{endpoint}*direction;
        errors(directionIndex)=norm(numeric-analytic)/max(norm(analytic),eps);
    end
    rows{endpoint}={w,fixedResidual,median(errors),max(errors)};
end
alignment=cell2table(vertcat(rows{:}),'VariableNames', ...
    {'w','fixedPointResidual','medianRelativeJvpError','maximumRelativeJvpError'});
writetable(alignment,fullfile(cfg.OutputRoot,'model_alignment_smoke.tsv'), ...
    'FileType','text','Delimiter','\t');
disp(alignment);
