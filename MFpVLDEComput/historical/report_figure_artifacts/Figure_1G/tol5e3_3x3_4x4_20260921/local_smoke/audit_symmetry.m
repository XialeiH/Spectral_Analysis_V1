function audit_symmetry(runRoot,outputRoot)
% Compare the existing representative reuse with evaluating every active row.
addpath(fullfile(runRoot,'repo','Utils'),'-begin');
addpath(fullfile(runRoot,'code'),'-begin');
addpath(fullfile(runRoot,'fast'),'-begin');
addpath(fullfile(outputRoot,'code'),'-begin');
maxNumCompThreads(12);
for fieldHC = [3 4]
    B = load(fullfile(runRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)));
    fast = h96_nn_set_condition_symmetry(h96_nn_prepare_base(B.ctx,'single', ...
        fullfile(runRoot,'fast','h96_models_double.mat')),100,0);
    fullFast = fast;
    fullFast.symRepresentative = (1:fast.activePairCount)';
    fullFast.symGroup = (1:fast.activePairCount)';
    fullFast.symClassCount = fast.activePairCount;
    for repeat = [1 4]
        state = B.initialState;
        if repeat > 3
            rng(1729+repeat,'twister');
            state.S = max(0,state.S.*exp(0.08*randn(size(state.S))-0.5*0.08^2));
            state.C = max(0,state.C.*exp(0.08*randn(size(state.C))-0.5*0.08^2));
            state.I = max(0,state.I.*exp(0.08*randn(size(state.I))-0.5*0.08^2));
        end
        for iteration = 1:5
            S=state.S; C=state.C; I=state.I;
            E=.6923*S+.3077*C;
            scale=L6Convert(E,fast.EKp)./E;
            Suse=S.*scale; Cuse=C.*scale; Iuse=InhMulp(I,fast.IKp);
            Euse=.6923*Suse+.3077*Cuse;
            recurrentE=fast.AE*[Suse;Cuse]; recurrentI=fast.AI*Iuse;
            field=reshape(Euse,fast.fieldRows,fast.fieldCols);
            filtered=conv2(padarray(field,[1 1],'circular'),fast.kernel,'same');
            l6Input=filtered(2:end-1,2:end-1);
            if isempty(fast.l6PP)
                l6=L6Convert(l6Input,fast.l6Pars);
            else
                l6=ppval(fast.l6PP,l6Input);
            end
            l6=min(max(l6(:),3),fast.l6Max)./3;
            reused=h96_nn_all_population_responses_symmetry(recurrentE,recurrentI,l6,fast);
            reference=h96_nn_all_population_responses_symmetry(recurrentE,recurrentI,l6,fullFast);
            relative=norm(reused(:)-reference(:))/max(norm(reference(:)),eps);
            maxAbs=max(abs(reused-reference),[],'all');
            mae=mean(abs(reused-reference),'all');
            fprintf('FIELD=%d REPEAT=%d STEP=%d REL=%.9g MAX_ABS=%.9g MAE=%.9g\n', ...
                fieldHC,repeat,iteration,relative,maxAbs,mae);
            state=struct('S',fast.p*reused(:,1)+fast.oneMinusP*S, ...
                'C',fast.p*reused(:,2)+fast.oneMinusP*C, ...
                'I',fast.p*reused(:,3)+fast.oneMinusP*I);
        end
    end
end
end

