function dataPaiCo = dataStructurePaiCo(proxies,instrumental,target,doBinning)

% Add some
% info/comment what this function does and
% how the format of the input parameters have to be like
% here

%% read in proxy data
for i = 1:numel(proxies)
    thisProxy.name = proxies(i).dataSetName; %Not needed, but might help identifying records
    thisProxy.data = proxies(i).paleoData_values';
    thisProxy.times = proxies(i).age';
    if doBinning
        for j = 1:numel(proxies(i).age)
            if j == 1
                thisProxy.lower(j) = proxies(i).age(j) - (proxies(i).age(j+1) - proxies(i).age(j)) / 2;
                thisProxy.upper(j) = (proxies(i).age(j) + proxies(i).age(j+1)) / 2;
            elseif j == numel(proxies(i).age)
                thisProxy.lower(j) = (proxies(i).age(j-1) + proxies(i).age(j)) / 2;
                thisProxy.upper(j) = proxies(i).age(j) + (proxies(i).age(j) - proxies(i).age(j-1)) / 2;
            else
                thisProxy.lower(j) = (proxies(i).age(j-1) + proxies(i).age(j)) / 2;
                thisProxy.upper(j) = (proxies(i).age(j) + proxies(i).age(j+1)) / 2;
            end
        end
    end
    thisProxy.locations = [proxies(i).geo_latitude, proxies(i).geo_longitude];
    dataPaiCo.proxy{i} = thisProxy;
    clear thisProxy;
end


%% read in instrumental data
dataPaiCo.instrumental.data = instrumental.data';
dataPaiCo.instrumental.times = instrumental.times';
if isfield(instrumental,'noisestd')
    dataPaiCo.instrumental.noisestd = instrumental.noisestd;
end
if isfield(instrumental,'timesReal')
    dataPaiCo.instrumental.timesReal = instrumental.timesReal;
end
if isfield(target,'timesReal')
    dataPaiCo.target.timesReal = target.timesReal;
end


%% read in target
dataPaiCo.target = target;