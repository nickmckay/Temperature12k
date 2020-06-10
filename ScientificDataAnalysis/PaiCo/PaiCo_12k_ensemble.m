    %% Get help on PaiCo

% help rtstruct
% data = rtstruct()
% help paico
% paicoinfo = paico('info');
% paicoinfo.info
% paicoinfo.options

clear
addpath(genpath('..'));
addpath(genpath('~/Dropbox/ChristophTemp12k'));



%% Read in LiPD files from DB and save the time series

load('/Users/npm4/Dropbox/Temp12kSerialization/fullEnsembles/temp12kEnsembles.mat')

%% Things to choose
doNoise = 0; %Choose true of false. If true, set doBinnig to false!
if doNoise
    doBinning = 0;
    binWidth = 1;
    annualiseBins = 0;
else
    doBinning = 1; %Choose true or false
    binWidth = 100; %Choose width of bins; 2000/binWidth should be a natural number
    if ~(floor(2000/binWidth) == 2000/binWidth)
        msg = 'Choose a binWidth s.t. 2000/binWidth is a natural number!';
        error(msg)
    end
    annualiseBins = 1;
end
nLatBands = 6; %Choose number of latitudinal bands

%% Define proxy records to use

TSG_ind = find(strcmp('temp12kEnsemble',{TS.paleoData_inCompilation}')); %index the compilation


TSG = TS(TSG_ind);

%% Do some manual corrections...

%get the annual, summerOnly and winterOnly proxy records
TSG_annual_ind = find(strncmp('annual',{TSG.interpretation1_seasonalityGeneral}',6)); %index all Temp12k annual records
TSG_summer_ind = find(strncmp('summerOnly',{TSG.interpretation1_seasonalityGeneral}',10)); %index all Temp12k summerOnly records
TSG_winter_ind = find(strncmp('winterOnly',{TSG.interpretation1_seasonalityGeneral}',10)); %index all Temp12k winterOnly records
TSG_annual = TSG(TSG_annual_ind); %all Temp12k annual records
TSG_summer = TSG(TSG_summer_ind); %all Temp12k SummerOnly records
TSG_winter = TSG(TSG_winter_ind); %all Temp12k WinterOnly records

TSG_big = TSG(union(TSG_annual_ind,union(TSG_summer_ind,TSG_winter_ind)));

TSG_incl_ind = unique([TSG_annual_ind; TSG_summer_ind; TSG_winter_ind]); %index all Temp12k annual, summerOnly and winterOnly records
TSG_incl = TSG(TSG_incl_ind); %all Temp12k annual, summerOnly and winterOnly records

%% Choose target
%this file must have/represent the same lat bands as the latBands matrix in PaiCoByLat.m !!
load('~/Dropbox/ChristophTemp12k/ChristophMatlab/targetMedian.mat');
instrTarget = targetMedian.CPS;
%instrTarget = targetMedian.PCR;
%instrTarget = targetMedian.CCA;
%instrTarget = targetMedian.GraphEM;
%instrTarget = targetMedian.AM;
%instrTarget = targetMedian.DA;


%% Run PaiCo for the specified proxy records and latitudinal bands

%% All
nens = 500;

txa = tic;
allPaiCoBig = cell(nens,1);
for nei = 1:length(allPaiCoBig)
    nei
    allPaiCoBig{nei} = PaiCoByLatEnsemble(TSG_big, nLatBands, doBinning, binWidth, instrTarget, doNoise, annualiseBins); %! result is in CE
end
toc(txa)


save bigPaiCo2.mat allPaiCoBig

%% Extract ensembles
latBins = ({'60-90N', '30-60N', '0-30N', '0-30S', '30-60S', '60-90S'});
for lb = 1:6
    em = nan(120,nens);
for ei = 1:nens
em(:,ei) = allPaiCoBig{ei}(lb).signal; 
end
csvwrite(['paicoEnsemble500-' latBins{lb} '.csv'],[1950-(allPaiCoBig{ei}(lb).times*binWidth-binWidth/2)' em])

end


%% Plot Big calib
signal_lat_weight = [];
figure(1); clf;
for ei = 1:100
    PaiCoBig = allPaiCoBig{ei};
    
    for i = 1:6
        weight = (1/(pi*(4-i)/6-pi*(3-i)/6))*(sin(pi*(4-i)/6)-sin(pi*(3-i)/6));
        signal_lat_weight = [signal_lat_weight;PaiCoBig(i).signal*weight];
    end
    PaiCoBigMean = mean(signal_lat_weight);
    if ~annualiseBins
        figure(1)
        plot(PaiCoBig(1).times,PaiCoBig(1).signal)
        hold on
        for i = 2:numel(PaiCoBig)
            plot(PaiCoBig(i).times,PaiCoBig(i).signal)
        end
        plot(PaiCoBig(1).times,PaiCoBigMean,'LineWidth',4)
        hline = refline(0,0);
        hline.Color = 'k';
        legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
        hold off
    else
        figure(1)
        plot(PaiCoBig(1).times*binWidth-binWidth/2,PaiCoBig(1).signal)
        hold on
        for i = 2:numel(PaiCoBig)
            plot(PaiCoBig(i).times*binWidth-binWidth/2,PaiCoBig(i).signal)
        end
        plot(PaiCoBig(1).times*binWidth-binWidth/2,PaiCoBigMean,'LineWidth',4)
        hline = refline(0,0);
        hline.Color = 'k';
        legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    end
end

%% Annual

txa = tic;
PaiCoAnnual = PaiCoByLat(TSG_annual, nLatBands, doBinning, binWidth, instrTarget, doNoise, annualiseBins); %! result is in CE
toc(txa)
%% Plot annual calib
signal_lat_weight = [];
for i = 1:6
    weight = (1/(pi*(4-i)/6-pi*(3-i)/6))*(sin(pi*(4-i)/6)-sin(pi*(3-i)/6));
    signal_lat_weight = [signal_lat_weight;PaiCoAnnual(i).signal*weight];
end
PaiCoAnnualMean = mean(signal_lat_weight);
if ~annualiseBins
    figure(1)
    plot(PaiCoAnnual(1).times,PaiCoAnnual(1).signal)
    hold on
    for i = 2:numel(PaiCoAnnual)
        plot(PaiCoAnnual(i).times,PaiCoAnnual(i).signal)
    end
    plot(PaiCoAnnual(1).times,PaiCoAnnualMean,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
else
    figure(1)
    plot(PaiCoAnnual(1).times*binWidth-binWidth/2,PaiCoAnnual(1).signal)
    hold on
    for i = 2:numel(PaiCoAnnual)
        plot(PaiCoAnnual(i).times*binWidth-binWidth/2,PaiCoAnnual(i).signal)
    end
    plot(PaiCoAnnual(1).times*binWidth-binWidth/2,PaiCoAnnualMean,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
end
%% Plot annual uncalib
signal_lat_weight = [];
for i = 1:6
    weight = (1/(pi*(4-i)/6-pi*(3-i)/6))*(sin(pi*(4-i)/6)-sin(pi*(3-i)/6));
    signal_lat_weight = [signal_lat_weight;PaiCoAnnual(i).uncalibSignal*weight];
end
PaiCoAnnualMeanUncalib = mean(signal_lat_weight);
if ~annualiseBins
    figure(4)
    plot(PaiCoAnnual(1).times,PaiCoAnnual(1).uncalibSignal)
    hold on
    for i = 2:numel(PaiCoAnnual)
        plot(PaiCoAnnual(i).times,PaiCoAnnual(i).uncalibSignal)
    end
    plot(PaiCoAnnual(1).times,PaiCoAnnualMeanUncalib,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
else
    figure(4)
    plot(PaiCoAnnual(1).times*binWidth-binWidth/2,PaiCoAnnual(1).uncalibSignal)
    hold on
    for i = 2:numel(PaiCoAnnual)
        plot(PaiCoAnnual(i).times*binWidth-binWidth/2,PaiCoAnnual(i).uncalibSignal)
    end
    plot(PaiCoAnnual(1).times*binWidth-binWidth/2,PaiCoAnnualMeanUncalib,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
end

%% Summer
txs = tic;
PaiCoSummer = PaiCoByLat(TSG_summer, nLatBands, doBinning, binWidth, instrTarget, doNoise, annualiseBins); %! result is in CE
toc(txs)
%% Plot summer calib
signal_lat_weight = [];
for i = 1:6
    weight = (1/(pi*(4-i)/6-pi*(3-i)/6))*(sin(pi*(4-i)/6)-sin(pi*(3-i)/6));
    signal_lat_weight = [signal_lat_weight;PaiCoSummer(i).signal*weight];
end
PaiCoSummerMean = mean(signal_lat_weight);
if ~annualiseBins
    figure(2)
    plot(PaiCoSummer(1).times,PaiCoSummer(1).signal)
    hold on
    for i = 2:numel(PaiCoSummer)
        plot(PaiCoSummer(i).times,PaiCoSummer(i).signal)
    end
    plot(PaiCoSummer(1).times,PaiCoSummerMean,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
else
    figure(2)
    plot(PaiCoSummer(1).times*binWidth-binWidth/2,PaiCoSummer(1).signal)
    hold on
    for i = 2:numel(PaiCoSummer)
        plot(PaiCoSummer(i).times*binWidth-binWidth/2,PaiCoSummer(i).signal)
    end
    plot(PaiCoSummer(1).times*binWidth-binWidth/2,PaiCoSummerMean,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
end
%% Plot summer uncalib
signal_lat_weight = [];
for i = 1:6
    weight = (1/(pi*(4-i)/6-pi*(3-i)/6))*(sin(pi*(4-i)/6)-sin(pi*(3-i)/6));
    signal_lat_weight = [signal_lat_weight;PaiCoSummer(i).uncalibSignal*weight];
end
PaiCoSummerMeanUncalib = mean(signal_lat_weight);
if ~annualiseBins
    figure(5)
    plot(PaiCoSummer(1).times,PaiCoSummer(1).uncalibSignal)
    hold on
    for i = 2:numel(PaiCoSummer)
        plot(PaiCoSummer(i).times,PaiCoSummer(i).uncalibSignal)
    end
    plot(PaiCoSummer(1).times,PaiCoSummerMeanUncalib,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
else
    figure(5)
    plot(PaiCoSummer(1).times*binWidth-binWidth/2,PaiCoSummer(1).uncalibSignal)
    hold on
    for i = 2:numel(PaiCoSummer)
        plot(PaiCoSummer(i).times*binWidth-binWidth/2,PaiCoSummer(i).uncalibSignal)
    end
    plot(PaiCoSummer(1).times*binWidth-binWidth/2,PaiCoSummerMeanUncalib,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
end

%% Winter
txw = tic;
PaiCoWinter = PaiCoByLat(TSG_winter, nLatBands, doBinning, binWidth, instrTarget, doNoise, annualiseBins); %! result is in CE
if ~annualiseBins
    figure(3)
    plot(PaiCoWinter(1).times,PaiCoWinter(1).signal)
    hold on
    for i = 2:numel(PaiCoWinter)
        plot(PaiCoWinter(i).times,PaiCoWinter(i).signal)
    end
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S'}, 'Location', 'best')
    hold off
else
    figure(3)
    plot(PaiCoWinter(1).times*binWidth-binWidth/2,PaiCoWinter(1).signal)
    hold on
    for i = 2:numel(PaiCoWinter)
        plot(PaiCoWinter(i).times*binWidth-binWidth/2,PaiCoWinter(i).signal)
    end
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S'}, 'Location', 'best')
    hold off
end
toc(txw)

%% TEST Recalibrate to other target
%% Calibrate to target
%instrTarget = targetMedian.CPS;
%instrTarget = targetMedian.PCR;
%instrTarget = targetMedian.CCA;
%instrTarget = targetMedian.GraphEM;
instrTarget = targetMedian.AM;
%instrTarget = targetMedian.DA;

%% Reconstruction target
if ~doBinning
    binWidth = 1;
end
if ~annualiseBins
    instrumental_forrecalib.data = smooth(instrTarget(i,:),binWidth,'loess');
    instrumental_forrecalib.times = 1:1:2000;
    instrumental_forrecalib.times = instrumental_forrecalib.times'; %using CE
else
    instrumental_forrecalib.data = [];
    for j = 1:(2000/binWidth)
        instrumental_forrecalib.data(j) = mean(instrTarget(i,(1+(j-1)*(binWidth)):(j*binWidth)));
    end
    instrumental_forrecalib.data = instrumental_forrecalib.data';
    instrumental_forrecalib.timesReal = (binWidth/2):binWidth:(2000-binWidth/2);
    instrumental_forrecalib.timesReal = instrumental_forrecalib.timesReal';
    instrumental_forrecalib.times = (instrumental_forrecalib.timesReal+binWidth/2)/binWidth; %using CE
end

signal_lat_weight = [];
signal_lat_unweight = [];
for i = 1:6
    tocalibrate = PaiCoAnnual(i).uncalibSignal;
    tocalibratepartmean = mean(tocalibrate(ismember(PaiCoAnnual(i).times, instrumental_forrecalib.times)));
    tocalibratepartstd = std(tocalibrate(ismember(PaiCoAnnual(i).times, instrumental_forrecalib.times)));
    calibrated = (tocalibrate-tocalibratepartmean)*std(instrumental_forrecalib.data)/tocalibratepartstd + mean(instrumental_forrecalib.data);
    
    weight = (1/(pi*(4-i)/6-pi*(3-i)/6))*(sin(pi*(4-i)/6)-sin(pi*(3-i)/6));
    signal_lat_weight = [signal_lat_weight;calibrated*weight];
    signal_lat_unweight = [signal_lat_unweight;calibrated];
end

PaiCoAnnualMeanRecalib = mean(signal_lat_weight);
if ~annualiseBins
    figure(7)
    plot(PaiCoAnnual(1).times,signal_lat_unweight(1,:))
    hold on
    for i = 2:6
        plot(PaiCoAnnual(i).times,signal_lat_unweight(i,:))
    end
    plot(PaiCoAnnual(1).times,PaiCoAnnualMeanRecalib,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
else
    figure(7)
    plot(PaiCoAnnual(1).times*binWidth-binWidth/2,signal_lat_unweight(1,:))
    hold on
    for i = 2:6
        plot(PaiCoAnnual(i).times*binWidth-binWidth/2,signal_lat_unweight(i,:))
    end
    plot(PaiCoAnnual(1).times*binWidth-binWidth/2,PaiCoAnnualMeanRecalib,'LineWidth',4)
    hline = refline(0,0);
    hline.Color = 'k';
    legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S', 'Mean'}, 'Location', 'best')
    hold off
end
