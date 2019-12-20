function [A,counts] = selectProxies(cache, indices, useFilecache)

indices = sort(indices);

if (useFilecache)
    % Rewind cache to beginning and fetch first proxy
    varcache(2,cache.idInd);
    comparisons = varcache(3,cache.idInd);
else 
    comparisons = cache.comparisons{1};
end

counts = zeros(cache.maxIndex, 2);

cacheIndex = 1;
for k = 1:numel(indices)
    while (cacheIndex < indices(k))    
        cacheIndex = cacheIndex + 1;
        if (useFilecache)
            comparisons = varcache(3,cache.idInd);
        else
            comparisons = cache.comparisons{cacheIndex};                
        end
    end
    
    signs = sign(comparisons);
    ind = sub2ind(size(counts), abs(comparisons), (3-signs)/2);
    counts(ind) = counts(ind) + 1;        
end


% Return only patterns with nonzero counts
mask = sum(counts,2)>0;
counts = counts(mask,:);
A = cache.A(mask,:);