function cache = inferProxies(times, proxy, tracker, useFilecache)

tracker('PAICO>pairs');

% Prefix tree stores the patterns for the proxy sample time frames, i.e.,
% different sets of target time indices covered by all possible proxy
% samples
% Later it is converted to an array of pattern lists that can be used to
% recover the time indices of different patterns
paico_mex();

nrProxy = 0;
for i = 1:numel(proxy)
    nrProxy = nrProxy + size(proxy{i}.data,1);
end

% Indices of pairwise comparioson for each proxy, pointing to rows in A
% Negative index means the relation needs to be inverted.
if (useFilecache)
    cache.idInd = varcache(0);
else
    cache.comparisons = cell(nrProxy,1);
end

% Handle records with the same resolution as the target differently as they
% populate the prefix tree and slow down the inference process
maxInd = numel(times)*(numel(times)-1)/2;
nnzA = 2*maxInd;

proxyInd = 0;

for pi = 1:numel(proxy)
    
    data = proxy{pi}.data;   
    
    inTargetRes = false;
    if (isfield(proxy{pi},'intargetres'))
        inTargetRes = proxy{pi}.intargetres;
    end
    
    if (~isfield(proxy{pi},'lower'))
        % Proxy record doesn't contain time frames for samples. Assume it
        % is annual.
        inTargetRes = true;
    end
    
    timeind = cell(numel(times),1);
    map = nan(size(timeind));
    pos = 0;
    
    if (inTargetRes)
        % Handle annual times separately
        
        % First sort times ascending
        [proxyTimes, ind] = sort(proxy{pi}.times, 'ascend');
        data = data(:,ind);
        
        % Only include samples that are in target
        [mask, loc] = ismember(proxyTimes, times);
        data = data(:,mask);
        
        % Find what are the target indices for the data in these proxies
        targetInd = loc(loc > 0);
                
        [indMin, indMaj] = find(tril(true(numel(targetInd)),-1));
        indices = (targetInd(indMin)-1).*(targetInd(indMin)-2)/2 + targetInd(indMaj);
                
        % Separate proxies to their own cells
        for i = 1:size(data,1)
            proxyComp = int32(sign(data(i,indMaj)-data(i,indMin)).*indices);
            proxyComp = proxyComp(abs(proxyComp) > 0);
            if (useFilecache)
                varcache(1,cache.idInd, proxyComp);
            else
                proxyInd = proxyInd + 1;
                cache.comparisons{proxyInd} = proxyComp;
            end 
        end
    
    else            
        % Fractions of years are rounded to closest integer year. Higher
        % resolutions than the target can not be reasonably handled with 
        % this implementation.
        lower = round(proxy{pi}.lower);
        upper = round(proxy{pi}.upper);

        % Sort data so that times are ascending
        [lower,ind] = sort(lower,'ascend');
        upper = upper(ind);
        data = data(:,ind);


        for i = 1:size(data,2)
            % Lower bound in time is inclusive, upper bound is exclusive.
            times1 = find(times >= lower(i) & times < upper(i));
            if (~isempty(times1)), 
                pos = pos + 1;
                timeind{pos} = times1;
                map(pos) = i;
                continue; 
            end
        end
    
        comparisons = int32(zeros(size(data,1),pos*(pos-1)/2));
        compPos = 0;

        for i = 1:pos
            times1 = timeind{i};

            for j = (i+1):pos
                times2 = timeind{j};

                signs = sign(data(:,map(i))-data(:,map(j)));
                % Skip comparison if signs are zero or values are NaN's
                if (all(isnan(signs) | (signs==0)))
                    continue;
                end

                ind = paico_mex(times1,times2);

                if (ind < 0) 
                    % Create new pattern
                    maxInd = maxInd + 1;
                    ind = maxInd;
                    nnzA = nnzA + numel(times1) + numel(times2);

                    paico_mex(times1,ind);
                    paico_mex(times2,ind);           
                end
                
                compPos = compPos + 1;
                comparisons(:,compPos) = ind*signs;
            end
        end
    
        % Separate proxies to their own cells
        for i = 1:size(comparisons,1)
            proxyComp = comparisons(i,1:compPos);
            proxyComp = proxyComp(abs(proxyComp) > 0);
            if (useFilecache)
                varcache(1,cache.idInd, proxyComp);
            else
                proxyInd = proxyInd + 1;
                cache.comparisons{proxyInd} = proxyComp;
            end        
        end

    end
    tracker('PAICO>pairs', pi, numel(proxy));
end

% Store measurements to conserve space when creating variables
cache.maxIndex = maxInd;
cache.nnzA = nnzA;
cache.nrTimes = numel(times);

% Convert prefix tree to arrays of patterns
paico_mex(-1);

% Create cache of coparisons matrix A
rows = zeros(nnzA,1);
cols = zeros(nnzA,1);
vals = zeros(nnzA,1);

% Add target resolved comparisons
[indMaj, indMin] = find(triu(true(numel(times)),1));
pos = numel(indMaj)*2;
rows(1:pos) = reshape([1;1]*(1:numel(indMaj)), pos, 1);
cols(1:pos) = reshape([indMaj';indMin'], pos, 1);
vals(1:pos) = reshape(([1;-1]./sqrt(2))*ones(1,numel(indMaj)), pos, 1);

% Add other resolved comparisons
for i = ((numel(times)-1)*numel(times)/2+1):maxInd
    [low high] = paico_mex(i);
    c = 1./sqrt(1/numel(low)+1/numel(high));
    
    ind = pos+(1:numel(low));
    rows(ind) = i;
    cols(ind) = low;
    vals(ind) = c/numel(low);
    pos = pos + numel(low);
    
    ind = pos+(1:numel(high));
    rows(ind) = i;
    cols(ind) = high;
    vals(ind) = -c/numel(high);
    pos = pos + numel(high);
end

cache.A = sparse(rows,cols,vals);
cache.nrProxy = nrProxy;

