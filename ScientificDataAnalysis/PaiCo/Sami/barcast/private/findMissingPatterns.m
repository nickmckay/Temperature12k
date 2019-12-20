function missingPatterns = findMissingPatterns(data, timePoints)
%% Recover the patterns of missing data
dataPatterns = cell(numel(data.proxy)+1,1);
timePatterns = zeros(timePoints, size(dataPatterns,1));

[dataPatterns{1} positions] = uniquePatterns(~isnan(data.instrumental.data'));
timePatterns(data.instrumental.timeIndices,1) = positions;

for j = 1:numel(data.proxy)
    nanMask = ~isnan(data.proxy{j}.data);
    if (nnz(nanMask) > 1)
        [patterns, positions] = uniquePatterns(nanMask');
        for i = 1:numel(patterns)
            % Indices to locations 
            locations = data.proxy{j}.locationIndices(patterns{i});
            [ind,~,pos] = unique(locations);
            pattern.dataIndices = patterns{i};
            pattern.locationIndices = ind;
            pattern.locationMap = double(repmat((1:numel(ind))',1,numel(pos)) == ...
                repmat(pos,numel(ind),1));
            patterns{i} = pattern;
        end
        dataPatterns{j+1} = patterns;
        timePatterns(data.proxy{j}.timeIndices,j+1) = positions;
    else
        % If no missing patterns in proxy, we can handle thigs faster
        [ind,~,pos] = unique(data.proxy{j}.locationIndices);
        pat = cell(1,1);
        pat{1}.locationIndices = ind;
        pat{1}.dataIndices = 1:size(data.proxy{j}.data,1);
        pat{1}.locationMap = double(repmat((1:numel(ind))',1,numel(pos)) == ...
            repmat(pos,numel(ind),1));
        dataPatterns{j+1} = pat;
        timePatterns(data.proxy{j}.timeIndices,j+1) = 1;
    end
end

[missingPatterns.timePatterns,~,missingPatterns.timeToPattern] = unique(timePatterns,'rows');
missingPatterns.dataPatterns = dataPatterns;

function [patterns positions] = uniquePatterns(data)
[rows,~,positions] = unique(data,'rows');
patterns = cell(size(rows,1),1);
for i = 1:size(rows,1)
    patterns{i} = find(rows(i,:));
end
