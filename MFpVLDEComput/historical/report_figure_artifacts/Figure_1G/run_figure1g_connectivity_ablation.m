function run_figure1g_connectivity_ablation(sourceRoot,outputRoot,fieldHC,repeatId)
% Compare exact sparse, precomputed-neighbor, and structure-agnostic dense products.
fieldHC=str2double(string(fieldHC)); repeatId=str2double(string(repeatId));
addpath(fullfile(sourceRoot,'fast'),'-begin'); addpath(fullfile(sourceRoot,'code'),'-begin');
threads=str2double(getenv('SLURM_CPUS_PER_TASK')); if isfinite(threads)&&threads>=1,maxNumCompThreads(threads);end
B=load(fullfile(sourceRoot,'bundles',sprintf('field_%02dHC.mat',fieldHC)),'ctx','initialState');
base=h96_nn_prepare_base(B.ctx,'single',fullfile(sourceRoot,'fast','h96_models_double.mat'));
matrices={base.AE,base.AI}; vectors={[B.initialState.S;B.initialState.C],B.initialState.I}; names=["Excitatory","Inhibitory"];
rows=cell(6,1); rowId=0;
for pathway=1:2
    A=matrices{pathway}; x=vectors{pathway}; reference=A*x;
    [neighborColumns,neighborWeights]=neighbor_form(A);
    neighborFn=@() sum(neighborWeights.*x(neighborColumns),2);
    sparseFn=@() A*x;
    [sparseMedian,sparseMean,sparseStd]=timing(sparseFn);
    [neighborMedian,neighborMean,neighborStd]=timing(neighborFn);
    neighborError=norm(neighborFn()-reference)/max(norm(reference),eps);
    rowId=rowId+1; rows{rowId}=make_row(fieldHC,repeatId,names(pathway),"Sparse", ...
        sparseMedian,sparseMean,sparseStd,0,0,nnz(A),size(A,1),size(A,2));
    rowId=rowId+1; rows{rowId}=make_row(fieldHC,repeatId,names(pathway),"LocalNeighbor", ...
        neighborMedian,neighborMean,neighborStd,0,neighborError,nnz(A),size(A,1),size(A,2));
    if fieldHC<=10
        timer=tic; denseA=full(A); conversion=toc(timer); denseFn=@() denseA*x;
        [denseMedian,denseMean,denseStd]=timing(denseFn); denseError=norm(denseFn()-reference)/max(norm(reference),eps);
    else
        conversion=NaN; denseMedian=NaN; denseMean=NaN; denseStd=NaN; denseError=NaN;
    end
    rowId=rowId+1; rows{rowId}=make_row(fieldHC,repeatId,names(pathway),"Dense", ...
        denseMedian,denseMean,denseStd,conversion,denseError,nnz(A),size(A,1),size(A,2));
end
T=vertcat(rows{:}); out=fullfile(outputRoot,'connectivity',sprintf('connectivity_%02dHC_repeat_%02d.tsv',fieldHC,repeatId));
writetable(T,out,'FileType','text','Delimiter','\t');
end

function [columns,weights]=neighbor_form(A)
[column,row,value]=find(A'); counts=accumarray(row,1,[size(A,1),1]); maxCount=max(counts);
starts=cumsum([0;counts(1:end-1)]); slot=(1:numel(value))'-repelem(starts,counts);
linear=sub2ind([size(A,1),maxCount],row,slot);
columns=ones(size(A,1),maxCount); weights=zeros(size(A,1),maxCount);
columns(linear)=column; weights(linear)=value;
end

function [med,avg,sd]=timing(fn)
fn(); samples=zeros(7,1); for k=1:numel(samples),timer=tic;fn();samples(k)=toc(timer);end
med=median(samples);avg=mean(samples);sd=std(samples);
end

function T=make_row(fieldHC,repeatId,pathway,representation,med,avg,sd,conversion,error,nonzeros,nRows,nColumns)
T=table(fieldHC,repeatId,pathway,representation,med,avg,sd,conversion,error,nonzeros,nRows,nColumns, ...
    'VariableNames',{'FieldHC','Repeat','Pathway','Representation','MedianMultiplySeconds', ...
    'MeanMultiplySeconds','StdMultiplySeconds','ConversionSeconds','RelativeError','Nonzeros','Rows','Columns'});
end
