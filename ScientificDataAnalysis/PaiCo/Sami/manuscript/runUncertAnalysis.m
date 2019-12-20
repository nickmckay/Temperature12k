function runUncertAnalysis

nrproxy = [10 20 50 100];
calibLength = 50;
autocorr = [.5 .9]';
snr = 0.4:0.4:3.2;
runs = 10;
targetLength = 500;
nrSamples = 10;

% Collect jobs to be used for parallel computing
dataCollection = cell(numel(nrproxy)*numel(snr)*runs,1);
pos = 0;
for ni = 1:numel(nrproxy)
    for si = 1:numel(snr)
        for ri = 1:runs
            t = zeros(numel(autocorr),targetLength);
            t(:,1) = randn(1,2)/2;
            for i = 2:targetLength
                t(:,i) = t(:,i-1).*autocorr + (1-autocorr).*randn(size(autocorr))/2;
            end
            t = sum(zscore(t,[],2),1)+randn;

            % Target data is stored only for printing purposes
            data.target.data = t;
            data.target.times = 1:targetLength;

            % Instrumental
            d = t((targetLength-calibLength+1):end);
            data.instrumental.data = d + randn(size(d))*std(d)/4; % SNR 3
            data.instrumental.noisestd = std(d)/4;
            data.instrumental.times = data.target.times((targetLength-calibLength+1):end);
            data.instrumental.locations = [1 1];


            % Proxy data
            proxyData = zeros(nrproxy(ni),targetLength);
            for i = 1:nrproxy(ni)
                d = t + randn(size(t))*std(t)/snr(si);
                proxyData(i,:) = d;
            end
            data.proxy{1}.data = proxyData;
            data.proxy{1}.locations = repmat([1 1],nrproxy(ni),1);
            data.proxy{1}.times = data.target.times;

            data.target.noisestd = std(t)/snr(si);

            pos = pos + 1;
            dataCollection{pos} = data;
        end
    end
end

paicoopts.regcov = 100;
paicoopts.maxIters = 1e2;
paicoopts.errorTolerance = 1e-8;
paicoopts.resample.number = nrSamples;
% paicoopts.resample.parallelize = true;

disp(['Total data: ' num2str(numel(dataCollection))]);
ppm = ParforProgMon('Uncertainty ', numel(dataCollection), 1, 300, 80);
tic;
% tracker('uncert');
index = fastrandperm(1:numel(dataCollection));
dataCollection = dataCollection(index);
uncertResult = cell(numel(dataCollection),1);
parfor i = 1:numel(dataCollection)
    data = dataCollection{i};
    result = [];
    
    % Store target data
    result.target = data.target.data;

    % PaiCo
    res = paico(data, paicoopts);
    result.PaiCo = res.signal;
    
    % Uncertainty coverage percentages
    result.uncert = mean(repmat(abs(result.PaiCo-result.target),nrSamples,1) >= abs(res.resample.signals-repmat(result.PaiCo,nrSamples,1)),1);
    
    uncertResult{i} = result;   
    
%     tracker('uncert',i,numel(dataCollection));
    ppm.increment();
end
[~,ind] = sort(index);
uncertResult = uncertResult(ind);

toc;
clear dataCollection;
ppm.delete()
mpath = fileparts(mfilename('fullpath'));
save([mpath filesep 'results' filesep 'uncertaintyExperiment_ols.mat']);


