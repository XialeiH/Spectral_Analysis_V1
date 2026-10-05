function compact = h96_nn_prepare_compact_config(page)
% Remove zero-weight pixel/LGN pairs before evaluating the h96 predictors.
compact = page;
[rowIndex, categoryIndex, pairWeight] = find(page.PixLGNCtgr);
linearIndex = rowIndex+(categoryIndex-1)*page.n;
pairCount = numel(rowIndex);

compact.mixRows = rowIndex;
compact.mixCategories = categoryIndex;
compact.mixAggregate = sparse(rowIndex, 1:pairCount, pairWeight, page.n, pairCount);
compact.staticCompactPages = zeros(pairCount, 5, 3, page.nnPrecision);
for population = 1:3
    compact.staticCompactPages(:,:,population) = ...
        page.staticPages(linearIndex,:,population);
end
compact.activePairCount = pairCount;
end
