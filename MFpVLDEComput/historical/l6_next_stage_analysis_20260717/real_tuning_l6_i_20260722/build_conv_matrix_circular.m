function K = build_conv_matrix_circular(kernel, Ny, Nx)
% Sparse circular convolution matrix matching the canonical conv2 path.

NyInput = double(Ny);
NxInput = double(Nx);
Ny = round(NyInput);
Nx = round(NxInput);
if ~isfinite(NyInput) || ~isfinite(NxInput) || Ny<1 || Nx<1 || ...
        abs(NyInput-Ny)>1e-8 || abs(NxInput-Nx)>1e-8
    error('RealTuning:MapSize','Ny and Nx must be positive integers.');
end

persistent cache
key = sprintf('k%dx%d_%dx%d',Ny,Nx,size(kernel,1),size(kernel,2));
if ~isempty(cache) && isfield(cache,key)
    K = cache.(key);
    return
end

n = Ny*Nx;
indexGrid = reshape(1:n,Ny,Nx);
[kernelHeight,kernelWidth] = size(kernel);
centerRow = ceil(kernelHeight/2);
centerColumn = ceil(kernelWidth/2);
rows = [];
columns = [];
values = [];
for row = 1:kernelHeight
    for column = 1:kernelWidth
        weight = kernel(row,column);
        if weight==0; continue; end
        source = circshift(indexGrid,[row-centerRow,column-centerColumn]);
        rows = [rows;indexGrid(:)]; %#ok<AGROW>
        columns = [columns;source(:)]; %#ok<AGROW>
        values = [values;weight*ones(n,1)]; %#ok<AGROW>
    end
end
K = sparse(rows,columns,values,n,n);
if isempty(cache); cache=struct(); end
cache.(key)=K;
end
