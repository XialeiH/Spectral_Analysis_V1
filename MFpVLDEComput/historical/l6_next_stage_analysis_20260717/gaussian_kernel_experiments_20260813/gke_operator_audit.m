function audit=gke_operator_audit(baseline,controlled,mapSize)
% Quantify what each connectivity intervention preserves and destroys.

names=["C_SS";"C_SC";"C_SI";"C_CS";"C_CC"; ...
    "C_CI";"C_IS";"C_IC";"C_II";"L6"];
rowSumRelativeError=zeros(numel(names),1);
meanRowSumRatio=zeros(numel(names),1);
controlledRowSumCV=zeros(numel(names),1);
translationError=zeros(numel(names),1);
rotation90Error=zeros(numel(names),1);
baselineRMSRadius=zeros(numel(names),1);
controlledRMSRadius=zeros(numel(names),1);
relativeFrobeniusChange=zeros(numel(names),1);
for index=1:numel(names)
    original=sparse(baseline.(names(index)));
    current=sparse(controlled.(names(index)));
    originalSums=full(sum(original,2));
    currentSums=full(sum(current,2));
    rowSumRelativeError(index)=norm(currentSums-originalSums)/max(norm(originalSums),eps);
    meanRowSumRatio(index)=mean(currentSums)/max(mean(originalSums),eps);
    controlledRowSumCV(index)=std(currentSums)/max(abs(mean(currentSums)),eps);
    relativeFrobeniusChange(index)=norm(current-original,'fro')/max(norm(original,'fro'),eps);
    translationError(index)=local_symmetry_error(current,mapSize,"translation");
    rotation90Error(index)=local_symmetry_error(current,mapSize,"rotation90");
    baselineRMSRadius(index)=local_rms_radius(original,mapSize);
    controlledRMSRadius(index)=local_rms_radius(current,mapSize);
end
audit=table(names,rowSumRelativeError,meanRowSumRatio,controlledRowSumCV, ...
    translationError,rotation90Error,baselineRMSRadius,controlledRMSRadius, ...
    relativeFrobeniusChange, ...
    'VariableNames',{'operator','rowSumRelativeError','meanRowSumRatio', ...
    'controlledRowSumCV','translationError','rotation90Error', ...
    'baselineRMSRadius','controlledRMSRadius','relativeFrobeniusChange'});
end

function value=local_symmetry_error(matrix,mapSize,type)
n=size(matrix,1);
grid=reshape(1:n,mapSize);
switch type
    case "translation"
        order=circshift(grid,[1 0]);
    case "rotation90"
        order=rot90(grid,1);
end
order=order(:);
value=norm(matrix(order,order)-matrix,'fro')/max(norm(matrix,'fro'),eps);
end

function rmsRadius=local_rms_radius(matrix,mapSize)
[rows,columns,values]=find(matrix);
[rowY,rowX]=ind2sub(mapSize,rows);
[columnY,columnX]=ind2sub(mapSize,columns);
dx=mod(columnX-rowX+floor(mapSize(2)/2),mapSize(2))-floor(mapSize(2)/2);
dy=mod(columnY-rowY+floor(mapSize(1)/2),mapSize(1))-floor(mapSize(1)/2);
weights=abs(values);
rmsRadius=sqrt(sum(weights.*(double(dx).^2+double(dy).^2))/sum(weights));
end
