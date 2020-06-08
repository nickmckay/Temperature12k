function [PaiCoReconByLat, PaiCoData] = PaiCoByLatEnsemble(proxyRec, nLatBands, doBinning, binWidth, instrTarget, doNoise, annualiseBins)

%% Check proxy interpretation direction and change values if necessary
for i = 1:numel(proxyRec)
    if strcmp(proxyRec(i).interpretation1_direction,'negative') %checking the climate interpretation direction
        proxyRec(i).paleoData_values = -1*proxyRec(i).paleoData_values; %flipping if interpretation direction negative
    end
end

%% Binning the records
% %OLD
% if doBinning
%     bins = 12050:-binWidth:-50; %from 12050 to -50 yr BP, i.e., from -10100 to 2000 CE
%     binsMidyear = bins(1:end-1) - binWidth/2;
%     for i = 1:numel(proxyRec)
%         values = [];
%         age = [];
%         for j = 1:(numel(bins)-1)
%             ageInd = find(bins(j) >= proxyRec(i).age & proxyRec(i).age > bins(j+1));
%             if numel(ageInd) >= 1
%                 values = [values;nanmean(proxyRec(i).paleoData_values(ageInd))];
%                 age = [age;binsMidyear(j)];
%             end
%         end
%         proxyRec(i).age = age;
%         proxyRec(i).paleoData_values = values;
%     end
% end

if doNoise
    for i = 1:numel(proxyRec)
        values = [];
        age = [];
        proxyRec(i).paleoData_valuesInfilled = [];
        proxyRec(i).ageInfilled = [];
        proxyRec(i).ageInfilledNoisy = [];
        proxyRec(i).ageOrig = proxyRec(i).age; %just for test purposes
        % for each year, choose the neares available value. Take the mean
        % in case 2 values are equally close
        proxyStart = min([12050,max(proxyRec(i).age)]);
        proxyEnd = max([-50,min(proxyRec(i).age)]);
        for k = proxyStart:-1:proxyEnd
            proxyRec(i).ageInfilled(proxyStart+1-k) = k;
            proxyRec(i).paleoData_valuesInfilled(proxyStart+1-k) = nanmean(proxyRec(i).paleoData_values(find(abs(proxyRec(i).age-k) == min(abs(proxyRec(i).age-k)))));
        end
        % now add noise
        proxyVar = nanvar(proxyRec(i).paleoData_values);
        wNoise = sqrt(proxyVar/20)*randn(size(proxyRec(i).paleoData_valuesInfilled)); %noise with with 5% of the variance fo the proxy record
        proxyRec(i).paleoData_valuesInfilledNoisy = proxyRec(i).paleoData_valuesInfilled + wNoise;
        
        proxyRec(i).age = proxyStart:-1:proxyEnd;
        proxyRec(i).age = proxyRec(i).age';
        proxyRec(i).paleoData_values = proxyRec(i).paleoData_valuesInfilledNoisy';
    end
end

% Annualise ("nearest neighbour") then bin to avoid large data gaps for low
% res records
if doBinning
    if doNoise
        msg = 'Set doBinning to false if doNoise is true!';
        error(msg)
    end
    
    
    
    bins = 12050:-binWidth:-50; %from 12050 to -50 yr BP, i.e., from -10100 to 2000 CE
    binsMidyear = bins(1:end-1) - binWidth/2;
    for i = 1:numel(proxyRec)
        
        
        %sample 1 member from the paleo ensemble
        npem = size(proxyRec(i).paleoData_values,2);
        if npem > 1
            proxyRec(i).paleoData_values = proxyRec(i).paleoData_values(:,randi(npem,1));
        end
        
        %sample 1 member from the age ensemble
        if isfield(proxyRec(i),'ageEnsemble')
            naem = size(proxyRec(i).ageEnsemble,2);
            if naem > 1
                proxyRec(i).age = proxyRec(i).ageEnsemble(:,randi(naem,1));
            end
            
        end
        
        
        
        values = [];
        age = [];
        proxyRec(i).paleoData_valuesInfilled = [];
        proxyRec(i).ageInfilled = [];
        proxyRec(i).ageOrig = proxyRec(i).age; %just for test purposes
        % for each year, choose the neares available value. Take the mean
        % in case 2 values are equally close
        proxyStart = min([12050,max(proxyRec(i).age)]);
        proxyEnd = max([-50,min(proxyRec(i).age)]);
        for k = proxyStart:-1:proxyEnd
            proxyRec(i).ageInfilled(proxyStart+1-k) = k;
            proxyRec(i).paleoData_valuesInfilled(proxyStart+1-k) = nanmean(proxyRec(i).paleoData_values(find(abs(proxyRec(i).age-k) == min(abs(proxyRec(i).age-k)))));
        end
        %do the binning. Calculate for each bin the mean of the values that
        %go into it and store the mid-year
        for j = 1:(numel(bins)-1)
            ageInd = find(bins(j) >= proxyRec(i).ageInfilled & proxyRec(i).ageInfilled > bins(j+1));
            if numel(ageInd) >= 1
                values = [values;nanmean(proxyRec(i).paleoData_valuesInfilled(ageInd))];
                age = [age;binsMidyear(j)];
            end
        end
        proxyRec(i).age = age;
        proxyRec(i).age = proxyRec(i).age*(-1)+1950; %convert to CE
        proxyRec(i).paleoData_values = values;
        if annualiseBins
            proxyRec(i).age = (proxyRec(i).age+binWidth/2)/binWidth;
        end
    end
end

%% Split analysis into latitudinal bands
latBandsWidth = 180 / nLatBands;
latBands = zeros(nLatBands,2); %Preallocate. Rows: bands; Column 1: min. lat; Column 2: max lat
for i = 1:nLatBands
    latBands(i,1) = (90 - (i-1)*latBandsWidth);
    latBands(i,2) = (90 - i*latBandsWidth);
end

for i = 1:nLatBands
    disp(i);
    proxiesInLatBand_ind = find(vertcat(proxyRec.geo_latitude) < latBands(i,1) & vertcat(proxyRec.geo_latitude) >= latBands(i,2));
    if numel(proxiesInLatBand_ind) == 0
        continue
    end
    proxiesInLatBand = proxyRec(proxiesInLatBand_ind);
    
    %% Reconstruction target
    if ~doBinning
        binWidth = 1;
    end
    if ~annualiseBins
        instrumental.data = smooth(instrTarget(i,:),binWidth,'loess');
        instrumental.times = 1:1:2000;
        instrumental.times = instrumental.times'; %using CE
        %instrumental.noisestd = .2224;
    else
        instrumental.data = [];
        for j = 1:(2000/binWidth)
            instrumental.data(j) = mean(instrTarget(i,(1+(j-1)*(binWidth)):(j*binWidth)));
        end
        instrumental.data = instrumental.data';
        instrumental.timesReal = (binWidth/2):binWidth:(2000-binWidth/2);
        instrumental.timesReal = instrumental.timesReal';
        instrumental.times = (instrumental.timesReal+binWidth/2)/binWidth; %using CE
    end
    
    %% Define target output
    if ~annualiseBins
        target.times = -10100:1:2000; %using CE
    else
        target.timesReal = (-10100+binWidth/2):binWidth:(2000-binWidth/2);
        target.times =(target.timesReal+binWidth/2)/binWidth;
    end
    
    %% Convert data into the data structure required for PaiCo
    dataPaiCoStucture = dataStructurePaiCo(proxiesInLatBand ,instrumental, target, doBinning);
    PaiCoData(i) = dataPaiCoStucture;
    
    %% Run PaiCo
    options.resample.number = 000;
    options.resample.parallelize = false;
    options.filecache = false;
    options.damping = [];
    options.errorTolerance = 1e-8;
    options.maxIters = 1e2;
    options.heuristicstart = true;
    options.regcov = 100;
    tracker();
    PaiCoReconByLat(i) = paico(dataPaiCoStucture, options, @tracker);
end

end