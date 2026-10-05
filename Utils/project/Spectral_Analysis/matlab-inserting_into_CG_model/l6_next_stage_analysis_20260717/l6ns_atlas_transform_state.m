function transformed = l6ns_atlas_transform_state(state,mapSize,rowShift,columnShift,reflectRows,reflectColumns)
% Apply one population-preserving spatial permutation to S, C, and I.

state=state(:);
n=prod(mapSize);
transformed=zeros(size(state),'like',state);
for population=1:3
    indices=(population-1)*n+(1:n);
    map=reshape(state(indices),mapSize);
    if reflectRows; map=flipud(map); end
    if reflectColumns; map=fliplr(map); end
    map=circshift(map,[rowShift columnShift]);
    transformed(indices)=map(:);
end
end
