function gke_smoke(setupFile,outputFile)
% Validate control invariants and baseline analytic-Jacobian parity.

loaded=load(setupFile,'setup');
setup=loaded.setup;
context=setup.Context;
catalog=gke_condition_catalog();
fprintf('LocalResponse: %s\n',which('LocalResponse_6D_MLP_prefAngle'));
fprintf('Iteration helper: %s\n',which( ...
    'LDEIteration_135FuncMain_CombDom_RealLGNL6_MLP6D_prefAngle'));
rows=cell(height(catalog),1);
for index=1:height(catalog)
    [operators,metadata,audit]=gke_build_operators(context,catalog.name(index));
    response=l6ns_controlled_phi(context.FixedPoint(:),context,operators);
    assert(all(isfinite(response)),'%s produced nonfinite values.',catalog.name(index));
    if metadata.ExpectedRowSumPreservation
        assert(metadata.MaximumL4RowSumRelativeError<1e-10, ...
            '%s did not preserve L4 row sums.',catalog.name(index));
    else
        assert(abs(metadata.MeanL4RowSumRatio-1)<1e-10, ...
            '%s did not preserve mean L4 strength.',catalog.name(index));
    end
    rows{index}=table(catalog.id(index),catalog.name(index), ...
        metadata.MaximumL4RowSumRelativeError,metadata.MeanL4RowSumRatio, ...
        metadata.MeanL4RowSumCV,metadata.MeanL4TranslationError, ...
        metadata.MeanL4Rotation90Error, ...
        'VariableNames',{'id','name','maximumRowSumRelativeError', ...
        'meanRowSumRatio','meanRowSumCV','translationError','rotation90Error'});
    fprintf('%02d %-38s row %.3e mean %.6f CV %.4f trans %.4f rot %.4f\n', ...
        index,catalog.name(index),metadata.MaximumL4RowSumRelativeError, ...
        metadata.MeanL4RowSumRatio,metadata.MeanL4RowSumCV, ...
        metadata.MeanL4TranslationError,metadata.MeanL4Rotation90Error);
    if catalog.name(index)=="swap_EI_ranges"
        assert(audit.controlledRMSRadius(1)<audit.baselineRMSRadius(1));
        assert(audit.controlledRMSRadius(3)>audit.baselineRMSRadius(3));
    end
end

[baselineOperators,~,~]=gke_build_operators(context,"baseline");
J=l6ns_controlled_jacobian(setup,baselineOperators,context.FixedPoint(:));
jacobianRelativeError=norm(J-setup.Pathway.JBaseline,'fro')/ ...
    max(norm(setup.Pathway.JBaseline,'fro'),eps);
assert(jacobianRelativeError<1e-12,'Baseline Jacobian parity failed.');
smokeSummary=vertcat(rows{:});
smokeSummary.baselineJacobianRelativeError(:)=jacobianRelativeError;
writetable(smokeSummary,outputFile,'FileType','text','Delimiter','\t');
fprintf('Smoke passed; baseline Jacobian relative error %.3e.\n',jacobianRelativeError);
end
