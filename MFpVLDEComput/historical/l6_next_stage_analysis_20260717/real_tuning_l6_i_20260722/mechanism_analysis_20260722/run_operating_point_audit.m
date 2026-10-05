function run_operating_point_audit()
% Separate pathway scaling, operating-point movement, and NN extension use.

paths = mechanism_initialize();
loaded = load(paths.SetupFile,'setup');
setup = loaded.setup;
context = setup.Context;
baseline = context.FixedPoint(:);
allCases = readtable(fullfile(paths.ResultsRoot,'all_case_summary.tsv'), ...
    'FileType','text','Delimiter','\t','TextType','string');

fprintf('Resolved predictor: %s\n',which('predict_pref6D_S'));
fprintf('Resolved extension: %s\n',which('h96_pref6D_domain_extended_response'));
fprintf('Resolved local response: %s\n',which('LocalResponse_6D_MLP_prefAngle'));

baselinePhiExact = l6ns_phi(baseline,0,context,[1 1],1,'true');
baselinePhiDiagnostic = mechanism_phi_variant(baseline,context,1,1,'extended');
phiMatchBaseline = norm(baselinePhiExact-baselinePhiDiagnostic)/ ...
    max(norm(baselinePhiExact),eps);
fprintf('Diagnostic Phi match at baseline: %.3e\n',phiMatchBaseline);
if phiMatchBaseline>1e-12
    error('Mechanism:PhiMismatch','Diagnostic Phi does not match the experiment map.');
end

specPathway = [repmat("L6",8,1);repmat("Inhibition",7,1)];
specBeta = [0;.15;.175;.20;.25;.275;.30;.40;0;-.05;-.075;-.10;-.125;-.15;-.20];
domainRows = cell(numel(specBeta)*2,1);
comparisonRows = cell(numel(specBeta),1);
for index = 1:numel(specBeta)
    pathway = specPathway(index);
    beta = specBeta(index);
    caseRow = local_case(allCases,pathway,beta);
    moved = local_load_state(paths.ResultsRoot,caseRow.taskId);
    gain6 = 1;
    gainI = 1;
    if pathway=="L6"; gain6 = 1+beta; else; gainI = 1+beta; end

    phiExact = l6ns_phi(moved,1-gain6,context,[1 1],gainI,'true');
    phiDiagnostic = mechanism_phi_variant(moved,context,gain6,gainI,'extended');
    phiError = norm(phiExact-phiDiagnostic)/max(norm(phiExact),eps);
    if phiError>1e-12
        error('Mechanism:PhiMismatch','Phi mismatch at %s beta %+.3f: %.3e.', ...
            pathway,beta,phiError);
    end
    domainRows{2*index-1} = mechanism_domain_stats( ...
        baseline,context,gain6,gainI,pathway,beta,'baseline state');
    domainRows{2*index} = mechanism_domain_stats( ...
        moved,context,gain6,gainI,pathway,beta,'moved equilibrium');

    if pathway=="L6"
        jFpp = setup.Pathway.JRest+gain6*setup.Pathway.J6+setup.Pathway.JI;
    else
        jFpp = setup.Pathway.JRest+setup.Pathway.J6+gainI*setup.Pathway.JI;
    end
    jDirectAtBaseline = real_tuning_true_jacobian( ...
        baseline,context,gain6,gainI);
    jBaselineMapAtMoved = real_tuning_true_jacobian(moved,context,1,1);
    jTrueAtMoved = real_tuning_true_jacobian(moved,context,gain6,gainI);
    maxFpp = mechanism_max_real(jFpp);
    maxDirectAtBaseline = mechanism_max_real(jDirectAtBaseline);
    maxBaselineMapAtMoved = mechanism_max_real(jBaselineMapAtMoved);
    maxTrueAtMoved = mechanism_max_real(jTrueAtMoved);
    comparisonRows{index} = table(pathway,beta,gain6,gainI, ...
        caseRow.relativeFixedPointShift,caseRow.maximumRateHz,maxFpp, ...
        maxDirectAtBaseline,maxBaselineMapAtMoved,maxTrueAtMoved, ...
        norm(jDirectAtBaseline-jFpp,'fro')/max(norm(jFpp,'fro'),eps), ...
        norm(jTrueAtMoved-jDirectAtBaseline,'fro')/ ...
            max(norm(jDirectAtBaseline,'fro'),eps),phiError, ...
        'VariableNames',{'pathway','beta','gain6','gainI', ...
        'relativeFixedPointShift','maximumRateHz','maxRealFpp', ...
        'maxRealDirectAtBaseline','maxRealBaselineMapAtMoved', ...
        'maxRealTrueAtMoved','relativeDirectVsFppJacobianDifference', ...
        'relativeMovedVsDirectJacobianDifference','phiReproductionError'});
    fprintf(['%s beta=%+.3f: maxRe FPP/direct0/base@moved/true@moved = ' ...
        '%.6f / %.6f / %.6f / %.6f\n'],pathway,beta,maxFpp, ...
        maxDirectAtBaseline,maxBaselineMapAtMoved,maxTrueAtMoved);
end

domainTable = vertcat(domainRows{:});
comparisonTable = vertcat(comparisonRows{:});
writetable(domainTable,fullfile(paths.OutputRoot,'domain_occupancy.tsv'), ...
    'FileType','text','Delimiter','\t');
writetable(comparisonTable,fullfile(paths.OutputRoot, ...
    'jacobian_operating_point_decomposition.tsv'),'FileType','text','Delimiter','\t');
save(fullfile(paths.OutputRoot,'operating_point_audit.mat'), ...
    'domainTable','comparisonTable','phiMatchBaseline','-v7.3');
local_plot(comparisonTable,domainTable,paths.FigureRoot);
end

function row = local_case(allCases,pathway,beta)
mask = allCases.pathway==pathway & abs(allCases.beta-beta)<1e-10;
if nnz(mask)~=1
    error('Mechanism:CaseLookup','Expected one %s beta %+.3f case.',pathway,beta);
end
row = allCases(mask,:);
end

function state = local_load_state(resultsRoot,taskId)
match = dir(fullfile(resultsRoot,sprintf('%03d_*',taskId),'case_result.mat'));
if numel(match)~=1
    error('Mechanism:CaseFile','Expected one result for task %d.',taskId);
end
loaded = load(fullfile(match.folder,match.name),'moved');
state = loaded.moved(:);
end

function local_plot(comparison,domain,figureRoot)
pathways = ["L6","Inhibition"];
colors = lines(4);
figure('Color','w','Position',[100 100 1260 500]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for p = 1:2
    nexttile;
    rows = comparison(comparison.pathway==pathways(p),:);
    [~,order] = sort(rows.beta);
    rows = rows(order,:);
    hFpp = plot(rows.beta,rows.maxRealFpp,'-o','LineWidth',1.5,'Color',colors(1,:)); hold on
    hDirect = plot(rows.beta,rows.maxRealDirectAtBaseline,'-s','LineWidth',1.5,'Color',colors(2,:));
    hBaseMoved = plot(rows.beta,rows.maxRealBaselineMapAtMoved,'-d','LineWidth',1.5,'Color',colors(3,:));
    hTrue = plot(rows.beta,rows.maxRealTrueAtMoved,'-^','LineWidth',1.5,'Color',colors(4,:));
    hThreshold = yline(1,'k--'); grid on
    xlabel('\beta'); ylabel('max Re(\lambda) of D\Phi');
    title(sprintf('%s: direct tuning versus operating-point movement',pathways(p)));
    if p==1
        legend([hFpp,hDirect,hBaseMoved,hTrue,hThreshold], ...
            {'FPP linear scaling','true gain at baseline state', ...
            'baseline map at moved state','true gain at moved equilibrium', ...
            'instability threshold'},'Location','northwest','FontSize',8);
    end
end
title(layout,'Why real tuning differs from fixed-point-preserving pathway scaling');
exportgraphics(gcf,fullfile(figureRoot,'01_operating_point_jacobian_decomposition.pdf'), ...
    'ContentType','vector');

figure('Color','w','Position',[100 100 1260 500]);
layout = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for p = 1:2
    nexttile;
    rows = domain(domain.pathway==pathways(p) & ...
        domain.stateLabel=="moved equilibrium",:);
    [~,order] = sort(rows.beta);
    rows = rows(order,:);
    yyaxis left
    plot(rows.beta,rows.fractionOutsideTrainingBox,'-o','LineWidth',1.5); hold on
    plot(rows.beta,rows.fractionFullExtension,'-s','LineWidth',1.5);
    ylabel('fraction of S/C/I sites'); ylim([-0.02 1.02]);
    yyaxis right
    plot(rows.beta,rows.rawVsExtendedRelativeResponseDifference,'-^','LineWidth',1.5);
    ylabel('relative raw-versus-extended response difference');
    xlabel('\beta'); grid on
    title(sprintf('%s: local-response domain at moved equilibrium',pathways(p)));
end
legend({'outside trained box','full extension','raw versus extended'}, ...
    'Location','southoutside','Orientation','horizontal');
title(layout,'Use of the h96 domain-extension branch');
exportgraphics(gcf,fullfile(figureRoot,'02_domain_extension_occupancy.pdf'), ...
    'ContentType','vector');
end
