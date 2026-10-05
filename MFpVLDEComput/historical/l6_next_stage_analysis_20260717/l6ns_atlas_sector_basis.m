function Q = l6ns_atlas_sector_basis(mapSize,rowWave,columnWave)
% Orthonormal basis for row-10 and column-20 translation sectors.

if ~isequal(mapSize,[40 40])
    error('The current sector basis expects a 40-by-40 map.');
end
rowsPerOrbit=mapSize(1)/4;
columnsPerOrbit=mapSize(2)/2;
n=prod(mapSize);
dimension=3*rowsPerOrbit*columnsPerOrbit;
nonzerosPerColumn=8;
rowIndices=zeros(nonzerosPerColumn*dimension,1);
columnIndices=zeros(nonzerosPerColumn*dimension,1);
values=complex(zeros(nonzerosPerColumn*dimension,1));
rowLambda=exp(1i*pi*rowWave/2);
columnLambda=(-1)^columnWave;
entry=0;
basisColumn=0;
for population=1:3
    populationOffset=(population-1)*n;
    for baseColumn=1:columnsPerOrbit
        for baseRow=1:rowsPerOrbit
            basisColumn=basisColumn+1;
            for columnShiftIndex=0:1
                column=baseColumn+columnShiftIndex*columnsPerOrbit;
                for rowShiftIndex=0:3
                    row=baseRow+rowShiftIndex*rowsPerOrbit;
                    entry=entry+1;
                    rowIndices(entry)=populationOffset+sub2ind(mapSize,row,column);
                    columnIndices(entry)=basisColumn;
                    values(entry)=rowLambda^(-rowShiftIndex)* ...
                        columnLambda^(-columnShiftIndex)/sqrt(8);
                end
            end
        end
    end
end
Q=sparse(rowIndices,columnIndices,values,3*n,dimension);
end
