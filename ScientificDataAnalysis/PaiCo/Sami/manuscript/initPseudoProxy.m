function initPseudoProxy(datafile)


nrproxy = [10 20 50 100];
% nrproxy = [10 50];
calibLength = 50;
autocorr = [.5 .9]';
snr = 0.4:0.4:3.2;
% snr = [0.4 2];
runs = 100;
targetLength = 500;

% Collect jobs to be used for parallel computing
dataCollection = cell(numel(nrproxy)*numel(snr)*2*runs,1);
pos = 0;
for ni = 1:numel(nrproxy)
    for si = 1:numel(snr)
        for nl = 1:2
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
                data.instrumental.data = d + randn(size(d))*std(d)/3; % SNR 3
                data.instrumental.times = data.target.times((targetLength-calibLength+1):end);
                data.instrumental.locations = [1 1];


                % Proxy data
                proxyData = zeros(nrproxy(ni),targetLength);
                for i = 1:nrproxy(ni)
                    d = t + randn(size(t))*std(t)/snr(si);

                    if (nl == 1)
                        x = [0 .25 .75 1];
                        y = [0 sort(rand(1,2)) 1];
                        d = interp1(x,y,(d-min(d))/(max(d)-min(d)),'cubic');                        
                    end
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
end

save(datafile);