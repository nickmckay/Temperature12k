
%Building Figure 6, Holocene temperature database
%July 2019
%% load the things

clear, close all

cols = brewermap(12,'Paired');
load TS.mat
load latBands2k30deg.mat;


%% Screen the database

mn = find(strncmpi('degC',{TS.paleoData_units}',1)); %index all calibrated
ma = find(strncmp('Temp12k',{TS.paleoData_inCompilation}',7)); %index the compilation

mb = find(strncmp('annual',{TS.interpretation1_seasonalityGeneral}',7)); %index all annual records
mc = find(strncmp('summerOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all SummerOnly records
me = find(strncmp('winterOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all WinterOnly records

mf = [mb; mc; me]; %index of annual, summer and winter

mg = intersect(mf,ma); %intersect seasons with compilation for a USE IN index

zog = intersect(mn,mg); % intersect x with the USE IN index to find calibrated index

hog = setdiff(mg,zog); %set difference between USE IN and calibrate to find uncalibrated index

%% grabbing all latitudes

latitude={TS.geo_latitude};
lat = (cell2mat(latitude))';

%% preallocating the latitudinal bands

Matrix(6).index = find(lat>=-90 & lat<-60);
Matrix(6).title = '90°S to 60°S';

Matrix(5).index = find(lat>=-60 & lat<-30);
Matrix(5).title = '60°S to 30°S';

Matrix(4).index = find(lat>=-30 & lat<05);
Matrix(4).title = '30°S to 0°S';

Matrix(3).index = find(lat>=0 & lat<30);
Matrix(3).title = '0°N to 30°N';

Matrix(2).index = find(lat>=30 & lat<60);
Matrix(2).title = '30°N to 60°N';

Matrix(1).index = find(lat>=60 & lat<90);
Matrix(1).title = '60°N to 90°N';

%% latitude steps
latstep = 90:-30:-90;

%% preallocate Bins

binStep=500; %SET BINSTEP
binEdges=[0 12000]; %SET PERIOD OF ANALYSIS (yr BP)
binVec=min(binEdges):binStep:max(binEdges);
binMid=binVec(1:end-1)+binStep/2;
binVec = binVec';
binMid = binMid';

plotEdges=reshape([binVec(1:end-1) binVec(2:end)]',[],1);

%% set the equal areas for griddng
normStart = 0; %SET START OF NORMALIZATION PERIOD (yr BP)
normEnd = 12000; %SET END OF NORMALIZATION PERIOD (yr BP)
[latgrid, longrid, regbound] = equalGridDeg(4000, -90);

%% zscore function with NaN values for uncalibrated records
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

%% Loop though each latitude 30° band and calculate latitudnal st*ff. Include # of records contributing to average.

for i = 1:6; %step through each latitude
    
    % Start the playing in a band (index all rows in latitudinal range)
    band = [latstep(i+1), latstep(i)];
    banding = find(cell2mat(latitude)>band(1,1) & cell2mat(latitude)< band(1,2)); %% specify the latitude cuttoff
    banding = banding';
    
    % index of records within latitude band (intersecting the records we
    % want in the latitudinal intervals
    uncal = intersect(hog,banding);
    cal = intersect(zog,banding); %calibrated
    
    %=======D Bin the calibrated records into matrices
    
    %do the things for one
    for j=1:length(cal) %step though each record in lat band
        %saving the record names in the latitude bands
        Matrix(i).namesCal{j,1} = TS(cal(j)).dataSetName;
        
        calLatData(j,1) = cell2mat({TS(cal(j)).geo_meanLat});
        calLonData(j,1) = cell2mat({TS(cal(j)).geo_meanLon});
        
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(cal(j)).age,TS(cal(j)).paleoData_values,binVec);
        if strcmp(TS(cal(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        calMat(:,j) = bin_mean;
        
    end
    
    clear j
    clear bin_mean
    
    %=======D Bin the uncalibrated records into matrices
    
    % do the things for the other
    if ~isempty(uncal)
        for j=1:length(uncal) %step though each record in lat band
            %saving the record names in the latitude bands
            Matrix(i).namesUncal{j,1} = TS(uncal(j)).dataSetName;
            
            uncalLatData(j,1) = cell2mat({TS(uncal(j)).geo_meanLat});
            uncalLonData(j,1) = cell2mat({TS(uncal(j)).geo_meanLon});
            
            [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(uncal(j)).age,TS(uncal(j)).paleoData_values,binVec);
            if strcmp(TS(uncal(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
                bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
            end
            
            uncalMat(:,j) = bin_mean;
            
        end
    else
        uncalMat(1:length(BinTime),1) = NaN;
    end
    
    %=======D getting sample depths
    
    dw = ~isnan(calMat);
    Matrix(i).calSampleDepth = sum(dw,2);
    
    dd = ~isnan(uncalMat);
    Matrix(i).uncalSampleDepth = sum(dd,2);
    
    %=======D mean and variance of 2k temperature bands (which are scaled to instrumental)--first bin to same as Holocene
    
    [bin2k, bin2ktime] = bin_x(bands2k(i).timeBP, bands2k(i).scaled_composite, binVec);
    good=find(bin2ktime>500 & bin2ktime<1500);
    at = nanmean(bin2k(good)); %23:73 = 520-1520 YR BP!!!!
    
    %=======D saving st*ff
    
    Matrix(i).calRecords = (calMat);
    Matrix(i).uncalRecords = (uncalMat);
    Matrix(i).time = binMid;
    Matrix(i).plotEdges = plotEdges;
    
    %=======D Gridding and compositing workflow.
    
    oMean=nanmean(calMat(binMid>normStart & binMid<normEnd,:),1); %using the mean for some period to compute anomalies %Default is to normalize over the whole 12000 year interval.
    oMat=repmat(oMean,length(BinTime),1); % turn it into a matrix for subtraction
    Tanom=calMat-oMat; % subract out the average
    [calComp,gridGroups,gridGroupsLat,gridGroupsLon,gridMean] = gridMat(Tanom,calLatData,calLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
    
    mx = nanmean(calComp(good));
    Matrix(i).calMedian = calComp-mx+at;
    
    %uncalibrated records
    if size(uncalMat,2)>1 %~isempty(uncal)
        Tunom=zscor_xnan(uncalMat); % zscores for uncalibrated rcords
        [uncalComp,unGridGroups,unGridGroupsLat,unGridGroupsLon,unGridMean] = gridMat(Tunom,uncalLatData,uncalLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
        Matrix(i).uncalMedian = uncalComp;
    else
        Tunom=zscor_xnan(uncalMat);
        Matrix(i).uncalMedian = Tunom;
    end
    
    %=======D Bootstrap ensembles
    
    m = bootstrp(1000,@nanmedian,gridMean');
    
    Matrix(i).calEnsemble = m';
    
    p = bootstrp(1000,@nanmean,unGridMean');
    Matrix(i).uncalEnsemble = p';
    
    %
    %clear st*ff for next loop
    clear bin_mean
    clear calMat
    clear uncalMat
    clear band
    clear banding
    clear calLatData
    clear calLonData
    clear uncalLatData
    clear uncalLonData
    i
end

bands = Matrix;
%% calculate weight

w(1,1) = (sind(90)-sind(60))/2;
w(1,2) = (sind(60)-sind(30))/2;
w(1,3) = (sind(30)-sind(0))/2;
w(1,4) = (sind(30)-sind(0))/2;
w(1,5) = (sind(60)-sind(30))/2;
w(1,6) = (sind(90)-sind(60))/2;

w1 = repmat(w, 24,1);

%% Converting from structure to column while doing st*ff
for k=1:6
    calLats(:,k) = bands(k).calMedian;
    uncalLats(:,k) = bands(k).uncalMedian;
    
    allSampleDepth(:,k) = sum([bands(k).calSampleDepth,bands(k).uncalSampleDepth],2);
    calSampleDepth(:,k) = bands(k).calSampleDepth;
    uncalSampleDepth(:,k) = bands(k).uncalSampleDepth;
    
    standardDev(:,k) = nanstd(bands(k).calEnsemble,0,2);
    uncalStd(:,k) =     nanstd(bands(k).uncalEnsemble,0,2);
    
    twoK(:,k) = bands2k(k).scaled_composite;
    
    %multiplying by weights
    twox(:,k) = w1(1,k)*bands2k(k).scaled_composite;
    instx(:,k) = w1(1,k)*bands2k(k).instTemp;
    cx(:,k) = w1(1,k)*bands(k).calMedian;
end



calGlobal = sum(cx,2);
% calGlobal  = wmean(calLats,w1,2);
% calGlobal = nanmean(calLats,2);
uncalGlobal = nanmean(uncalLats,2);
calSampleGlobal = nansum(calSampleDepth,2);
uncalSampleGlobal = nansum(uncalSampleDepth,2);
allSampleGlobal = nansum(allSampleDepth,2);
allStdGlobal = nanmean(standardDev,2);
uncalStdGlobal = nanmean(uncalStd,2);
twoThousand = sum(twox,2);%nanmean(twoK,2);
istruMental = sum(instx,2);

%%

% bands = Matrix;
%
%
% Fig6_plot
%
% % save the j*nk
save bandsFig6.mat bands

%% Plot
style_l = {'FontName','Heveltica','FontSize',12,'FontWeight','bold'};
figure(1)
clf

subplot(2,1,1);

[ax,h1,h2] = plotyy(bands(1).time,calGlobal,bands(1).time,uncalGlobal,@line,@line); hold on
hold(ax(1))
hold(ax(2))
hold on
plot(ax(1),bands2k(1).timeBP,twoThousand)
set(h1,'color',rgb('brown'),'linewidth',2)
set(ax(1),'XColor' , [.3 .3 .3], 'YColor', rgb('brown'), 'LineWidth', 1);
set(ax(1),'YAxisLocation','left','TickDir','out','YMinorTick','off')
set(get(ax(1),'Ylabel'),'String','temperature (°C)')
set(ax(1),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
set(ax(1),'Xdir','reverse','Xlim',[-100 12000])

%set stuff for line plots
set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
set(ax(2),'Xdir','reverse','Xlim',[-100 12000])
set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(get(ax(2),'Ylabel'),'String','temperature (zscore)');

set(h2,'color',rgb('darkslategray'),'linewidth',2)
set(ax(2),'XColor' , [.3 .3 .3], 'YColor', rgb('darkslategray'), 'LineWidth', 1);
set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(ax(1),'XMinorTick','on','YMinorTick','on', 'YGrid','on')

%plotting confidence intervals
area_fill(bands(1).time',[calGlobal+allStdGlobal]',[calGlobal-allStdGlobal]',rgb('brown'),rgb('brown'),1,0.2);


filled=[[uncalGlobal+uncalStdGlobal]',fliplr([uncalGlobal-uncalStdGlobal]')];
xpoints=[bands(1).time',fliplr(bands(1).time')];

fillhandle=fill(xpoints,filled,rgb('darkslategray'),'Parent',ax(2));%plot the data

set(fillhandle,'EdgeColor',rgb('darkslategray'),'FaceAlpha',0.2,'EdgeAlpha',0.2);%set edge color


subplot(2,1,2);

%records and sambple bars
[ax,h1,h2] = plotyy(bands(1).time,allSampleGlobal,bands(1).plotEdges,reshape([calSampleGlobal calSampleGlobal]',[],1),@bar,@line); hold on
hold(ax(1))
hold(ax(2))
%other records
plot(ax(2),bands(1).plotEdges,reshape([uncalSampleGlobal uncalSampleGlobal]',[],1), 'color', rgb('darkslategray'), 'linewidth',2);
hold on

%set stuff for bar graph
set(h1,'facecolor',rgb('Gainsboro'),'edgecolor','none')
set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','left','TickDir','out','YMinorTick','off')
set(ax(1),'yLim',[0 650])
set(ax(2),'yLim',[0 650])
set(ax(2),'Ytick',200:200:600)
set(ax(1),'Ytick',200:200:600)
set(ax(1),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
set(ax(1),'Xdir','reverse','Xlim',[-100 12000])
set(ax(1),'XMinorTick','on','YMinorTick','on', 'YGrid','on')

%set stuff for line plots
set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
set(ax(2),'Xdir','reverse','Xlim',[-100 12000])
set(get(ax(2),'Ylabel'),'String','#records');
set(get(ax(1),'Ylabel'),'String','total #records');
set(get(ax(2),'Ylabel'),style_l{:}),
set(gcf,'CurrentAxes',ax(2))
set(h2,'color',rgb('brown'),'linewidth',2)
set(ax(2),'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')