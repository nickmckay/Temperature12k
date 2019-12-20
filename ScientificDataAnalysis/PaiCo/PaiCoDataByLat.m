function [OrigData, NewData, PaiCoData] = PaiCoDataByLat(proxyRec, nLatBands, doBinning, binWidth, instrTarget)

%% Check proxy interpretation direction and change values if necessary
for i = 1:numel(proxyRec)
    if strcmp(proxyRec(i).interpretation1_direction,'negative') %checking the climate interpretation direction
            proxyRec(i).paleoData_values = -1*proxyRec(i).paleoData_values; %flipping if interpretation direction negative
    end
end

%% Binning the records
% %OLD
% if doBinning
%     bins = 11950:-binWidth:-50; %from 11950 to -50 yr BP, i.e., from -10000 to 2000 CE
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

% Annualise ("nearest neighbour") then bin to avoid large data gaps for low
% res records
if doBinning
    bins = 11950:-binWidth:-50; %from 11950 to -50 yr BP, i.e., from -10000 to 2000 CE
    binsMidyear = bins(1:end-1) - binWidth/2;
    for i = 1:numel(proxyRec)
        values = [];
        age = [];
        proxyRec(i).paleoData_valuesInfilled = [];
        proxyRec(i).ageInfilled = [];
        proxyRec(i).ageOrig = proxyRec(i).age; %just for test purposes
        proxyRec(i).valuesOrig = proxyRec(i).paleoData_values; %just for test purposes
        % for each year, choose the neares available value. Take the mean
        % in case 2 values are equally close
        proxyStart = min([11950,max(proxyRec(i).age)]);
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
        proxyRec(i).paleoData_values = values;
        
        %% Save new and original data
        OrigData(i).age = proxyRec(i).ageOrig;
        OrigData(i).values = proxyRec(i).valuesOrig;
        NewData(i).age = proxyRec(i).age;
        NewData(i).values = proxyRec(i).paleoData_values;
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
    instrumental.data = smooth(instrTarget(i,:),binWidth,'loess');
    instrumental.times = 1:1:2000;
    instrumental.times = instrumental.times'; %using CE
    %instrumental.noisestd = .2224;
    
    %% Define target output
    target.times = -10000:1:2000; %using CE
    
    %% Convert data into the data structure required for PaiCo
    dataPaiCoStucture = dataStructurePaiCo(proxiesInLatBand ,instrumental, target);
    PaiCoData(i) = dataPaiCoStucture;

end

end