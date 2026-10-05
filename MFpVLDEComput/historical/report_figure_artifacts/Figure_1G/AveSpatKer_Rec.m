function Cmean = AveSpatKer_Rec(source,nHCIn,nHCOutX,nHCOutY,nPixX,nPixY,varargin)
% Sparse periodic equivalent of the legacy large-field connectivity builder.
% Figure 1G uses Parax=1 and Paray=0, so ocular redistribution is identity.
if isempty(varargin)
    paraX = 1; paraY = 0;
else
    paraX = varargin{1}; paraY = varargin{2};
end
assert(paraX==1 && paraY==0,'Sparse benchmark builder requires Parax=1, Paray=0.');
assert(nHCOutX==nHCOutY, ...
    'Sparse benchmark builder requires a square HC field.');
assert(size(source,1)==size(source,2));
assert(size(source,1)==nHCIn*nPixX*nHCIn*nPixY);

kernelSamples = zeros(2*nPixY,2*nPixX,size(source,1));
for pixel = 1:size(source,1)
    px = ceil(pixel/(nHCIn*nPixY));
    py = mod(pixel,nHCIn*nPixY);
    if py==0, py=nHCIn*nPixY; end
    connMap = reshape(source(pixel,:).',nHCIn*nPixY,nHCIn*nPixX);
    centerHC = floor((nHCIn-1)/2);
    center = connMap(centerHC*nPixY+(1:nPixY),centerHC*nPixX+(1:nPixX));
    horizontal = connMap(centerHC*nPixY+(1:nPixY),:);
    vertical = connMap(:,centerHC*nPixX+(1:nPixX));
    extended = zeros((nHCIn+2)*nPixY,(nHCIn+2)*nPixX);
    extended(nPixY+(1:nHCIn*nPixY),nPixX+(1:nHCIn*nPixX)) = connMap;
    extended(1:nPixY,1:nPixX) = center;
    extended((nHCIn+1)*nPixY+(1:nPixY),1:nPixX) = center;
    extended(1:nPixY,(nHCIn+1)*nPixX+(1:nPixX)) = center;
    extended((nHCIn+1)*nPixY+(1:nPixY),(nHCIn+1)*nPixX+(1:nPixX)) = center;
    extended(1:nPixY,nPixX+(1:nHCIn*nPixX)) = horizontal;
    extended((nHCIn+1)*nPixY+(1:nPixY),nPixX+(1:nHCIn*nPixX)) = horizontal;
    extended(nPixY+(1:nHCIn*nPixY),1:nPixX) = vertical;
    extended(nPixY+(1:nHCIn*nPixY),(nHCIn+1)*nPixX+(1:nPixX)) = vertical;
    pxNew = px+nPixX; pyNew = py+nPixY;
    kernelSamples(:,:,pixel) = extended(pyNew-nPixY:pyNew+nPixY-1, ...
        pxNew-nPixX:pxNew+nPixX-1);
end
kernel = mean(kernelSamples,3);
tmp = zeros(2*nPixY,2*nPixX,2);
tmp(2:end,2:end,1) = kernel(2:end,2:end);
tmp(2:end,2:end,2) = kernel(end:-1:2,2:end);
kernel = mean(tmp,3);
tmp(:) = 0;
tmp(2:end,2:end,1) = kernel(2:end,end:-1:2);
tmp(2:end,2:end,2) = kernel(2:end,2:end);
kernel = mean(tmp,3);
tmp4 = zeros(2*nPixY,2*nPixX,4);
tmp4(2:end,2:end,1) = kernel(2:end,2:end);
tmp4(2:end,2:end,2) = kernel(2:end,2:end).';
kernel = mean(tmp4(:,:,1:2),3);
tmp4(2:end,2:end,3) = rot90(kernel(2:end,2:end),2);
tmp4(2:end,2:end,4) = kernel(2:end,2:end);
kernel = mean(tmp4,3);

[ky,kx,value] = find(kernel);
support = numel(value);
height = nHCOutY*nPixY;
width = nHCOutX*nPixX;
n = height*width;
rows = repelem((1:n).',support,1);
cols = zeros(n*support,1);
values = repmat(value,n,1);
for pixel = 1:n
    px = ceil(pixel/height);
    py = mod(pixel,height);
    if py==0, py=height; end
    yRange = mod(py-nPixY:py+nPixY-1,height);
    yRange(yRange==0) = height;
    xRange = mod(px-nPixX:px+nPixX-1,width);
    xRange(xRange==0) = width;
    linearMap = yRange(:) + (xRange(:).'-1)*height;
    block = (pixel-1)*support+(1:support);
    cols(block) = linearMap(sub2ind(size(linearMap),ky,kx));
end
Cmean = sparse(rows,cols,values,n,n,n*support);
end
