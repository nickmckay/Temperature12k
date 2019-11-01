%% Building Figure 8, Temperature 12k. Global composite
%Cody Routson, July 2019
%% load the things

clear, close all

load TS.mat
load PAGES_multiMethodMeadian.txt;
load grid.mat;

%% Screen the database

mn = find(strncmpi('degC',{TS.paleoData_units}',1)); %index all calibrated
ma = find(strcmp('Temp12k',{TS.paleoData_inCompilation}')); %index the compilation

mb = find(strncmp('annual',{TS.interpretation1_seasonalityGeneral}',7)); %index all annual records
mc = find(strncmp('summerOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all SummerOnly records
me = find(strncmp('winterOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all WinterOnly records

mf = [mb; mc; me]; %index of annual, summerOnly and winterOnly

mg = intersect(mf,ma); %seasons and Temp12k index

gm = intersect(mb,ma); %annual and Temp12k index

zog = intersect(mn,mg); % finding calibrated index

dog = intersect(mn,gm); % finding calibrated annual index

hog = setdiff(mg,zog); % finding uncalibrated index

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

%% Interval over which to subtract the mean (calibrated records) and normalize (uncalibrated records)
normStart = 0; %SET START OF NORMALIZATION PERIOD (yr BP)
normEnd = 12000; %SET END OF NORMALIZATION PERIOD (yr BP)

%% equal area latitude and longitude grid
latgrid = grid.latgrid;
longrid = grid.longrid;
%% calculating and saving the grid
% [latgrid, longrid, regbound] = equalGridDeg(4000, -90);
% grid.latgrid = latgrid;
% grid.longrid = longrid;
% grid.grebound = regbound;
% save grid.mat grid;


%% zscore function with NaN values for uncalibrated records
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

%% Loop though each latitude 30° band and calculate latitudnal st*ff. Include # of records contributing to average.

for i = 1:6; %step through each 30° latitude band
    
    % Start the playing in a band (index all rows in latitudinal range)
    band = [latstep(i+1), latstep(i)];
    banding = find(cell2mat(latitude)>band(1,1) & cell2mat(latitude)< band(1,2)); %% specify the latitude cuttoff
    banding = banding';
    
    % index of records within latitude band (intersecting the records we
    % want in the latitudinal intervals
    uncal = intersect(hog,banding);
    cal = intersect(zog,banding); %calibrated
    calAnn = intersect(dog,banding); %calibrated
    
    %=======D Bin the calibrated records into matrices
  
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
    
    %=======D Bin the calibrated ANNUAL records into matrices
  
    for j=1:length(calAnn) %step though each record in lat band
        %saving the record names in the latitude bands
        Matrix(i).namesCalAnn{j,1} = TS(calAnn(j)).dataSetName;
        
        calAnnLatData(j,1) = cell2mat({TS(calAnn(j)).geo_meanLat});
        calAnnLonData(j,1) = cell2mat({TS(calAnn(j)).geo_meanLon});
        
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(calAnn(j)).age,TS(calAnn(j)).paleoData_values,binVec);
        if strcmp(TS(calAnn(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        calAnnMat(:,j) = bin_mean;
        
    end
    
    clear j
    clear bin_mean
    
    %=======D Bin the uncalibrated records into matrices
    
    if ~isempty(uncal)
        for j=1:length(uncal) %step though each record in lat band
            %saving the record names in the latitude bands
            Matrix(i).namesUncal{j,1} = TS(uncal(j)).dataSetName;
            
            %collecting location information
            uncalLatData(j,1) = cell2mat({TS(uncal(j)).geo_meanLat});
            uncalLonData(j,1) = cell2mat({TS(uncal(j)).geo_meanLon});
            
            %binning
            [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(uncal(j)).age,TS(uncal(j)).paleoData_values,binVec);
            if strcmp(TS(uncal(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
                bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
            end
            
            uncalMat(:,j) = bin_mean;
            
        end
    else
        uncalMat(1:length(BinTime),1) = NaN;
    end
    
   

    %=======D saving st*ff
    
    Matrix(i).calRecords = (calMat);
    Matrix(i).uncalRecords = (uncalMat);
    Matrix(i).time = binMid;
    Matrix(i).plotEdges = plotEdges;
    
    %=======D Gridding and compositing workflow.
    
    % making anomalies gridding calibrated records
    oMean=nanmean(calMat(binMid>normStart & binMid<normEnd,:),1); %using the mean for some period to subtract from each record %Default is to subtract the mean of the whole 12000 year interval.
    oMat=repmat(oMean,length(BinTime),1); % turn it into a matrix for subtraction
    Tanom=calMat-oMat; % subract out the average
    [calComp,gridGroups,gridGroupsLat,gridGroupsLon,gridMean] = gridMat(Tanom,calLatData,calLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
    
    good = find(binMid>=500 & binMid<=1500);
    
    mx = nanmean(calComp(good)); %Finding the mean of the composite (500-1500 BP)
    Matrix(i).calMedian = calComp-mx;  %subtract the Holocene Mean of the composite
    
    % gridding calibrated ANNUAL records
    
    oMean2=nanmean(calAnnMat(binMid>normStart & binMid<normEnd,:),1); %using the mean for some period to compute anomalies %Default is to normalize over the whole 12000 year interval.
    oMat2=repmat(oMean2,length(BinTime),1); % turn it into a matrix for subtraction
    Tanom2=calAnnMat-oMat2; % subract out the average
    [calAnnComp,gridGroups2,gridGroupsLat2,gridGroupsLon2,gridMean2] = gridMat(Tanom2,calAnnLatData,calAnnLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
    
    mx2 = nanmean(calAnnComp(good)); %Finding the mean for some interval of the composite (500-1500 BP)
    Matrix(i).calAnnMedian = calAnnComp-mx2;  %subtract that mean of the composite
    
    %gridding uncalibrated records if they exist
    if size(uncalMat,2)>1 %
        Tunom=zscor_xnan(uncalMat); % zscores for uncalibrated rcords
        [uncalComp,unGridGroups,unGridGroupsLat,unGridGroupsLon,unGridMean] = gridMat(Tunom,uncalLatData,uncalLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
        Matrix(i).uncalMedian = uncalComp;
    else
        Tunom=zscor_xnan(uncalMat);
        Matrix(i).uncalMedian = Tunom;
    end
    
     %=======D getting sample depths
    
    dw = ~isnan(Tanom);
    Matrix(i).calSampleDepth = sum(dw,2);
    
    dg = ~isnan(Tanom2);
    Matrix(i).calAnnSampleDepth = sum(dg,2);
    
    dd = ~isnan(Tunom);
    Matrix(i).uncalSampleDepth = sum(dd,2);
    
    %=======D Bootstrap ensembles
    
    %calibrated records
    m = bootstrp(1000,@nanmedian,gridMean');
    Matrix(i).calEnsemble = m';
    
    %calibrated annual
    m2 = bootstrp(1000,@nanmedian,gridMean2');
    Matrix(i).calAnnEnsemble = m2';
    
    %uncalibrted
    p = bootstrp(1000,@nanmean,unGridMean');
    Matrix(i).uncalEnsemble = p';
    
    
    %clear st*ff for next loop
    clear bin_mean
    clear calMat
    clear calAnnMat
    clear uncalMat
    clear band
    clear banding
    clear calLatData
    clear calLonData
    clear uncalLatData
    clear uncalLonData
    clear calAnnLatData
    clear calAnnLonData
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

%% Converting from structure to column while doing st*ff. Sample depths, weighting, yata yata
for k=1:6
    calLats(:,k) = bands(k).calMedian;
    uncalLats(:,k) = bands(k).uncalMedian;
    
    allSampleDepth(:,k) = sum([bands(k).calSampleDepth,bands(k).uncalSampleDepth],2);
    calSampleDepth(:,k) = bands(k).calSampleDepth;
    calAnnSampleDepth(:,k) = bands(k).calAnnSampleDepth;
    uncalSampleDepth(:,k) = bands(k).uncalSampleDepth;
    
    standardDev(:,k) = nanstd(bands(k).calEnsemble,0,2);
    standardDevAnn(:,k) = nanstd(bands(k).calAnnEnsemble,0,2);
    uncalStd(:,k) =     nanstd(bands(k).uncalEnsemble,0,2);
    
    %multiplying by weights
    cx(:,k) = w1(1,k)*bands(k).calMedian;
    cx2(:,k) = w1(1,k)*bands(k).calAnnMedian;
    cx3(:,k) = w1(1,k)*bands(k).uncalMedian;
end


calGlobal = sum(cx,2); %averaging (summing weighted latitudinal bands) the global calibrated composite
calAnnGlobal = sum(cx2,2); %averaging (summing weighted latitudinal bands) the global Annual calibrated composite
uncalGlobal = nansum(cx3,2); %averaging (summing weighted latitudinal bands) the global uncalibrated composite

%sample depths for the various things
calSampleGlobal = nansum(calSampleDepth,2);
calAnnSampleGlobal = nansum(calAnnSampleDepth,2);
uncalSampleGlobal = nansum(uncalSampleDepth,2);
allSampleGlobal = nansum(allSampleDepth,2);
allStdGlobal = nanmean(standardDev,2);
calAnnStdGlobal = nanmean(standardDevAnn,2);
uncalStdGlobal = nanmean(uncalStd,2);

%% subtracting the correct median from the pages data
p2k.time = PAGES_multiMethodMeadian(:,1);
p2k.tmp = PAGES_multiMethodMeadian(:,2);

gd = find(p2k.time>=500 & p2k.time<=1500);

gt = nanmean(p2k.tmp(gd));
p2k.tmp = p2k.tmp-gt;

%bin the pages
[bin_mean2k, BinTime2k, sem, bin_sum] =  bin_x(p2k.time,p2k.tmp,binVec);

%% save some j*nk
% save bandsFig6.mat bands

%% Plot
style_l = {'FontName','Heveltica','FontSize',12,'FontWeight','bold'};
figure(1)
clf

%timeseries subplot
subplot(2,1,1);
[ax,h1,h2] = plotyy(bands(1).time,calGlobal,bands(1).time,uncalGlobal,@line,@line); hold on
hold(ax(1))
hold(ax(2))
hold on
plot(ax(1),bands(1).time,calAnnGlobal, 'color', rgb('black'), 'linewidth',2);
hold on
plot(ax(1),BinTime2k,bin_mean2k, 'color', rgb('red'), 'linewidth',2);
hold on
plot(ax(1),p2k.time,p2k.tmp, 'color', rgb('red'), 'linewidth',1);

hold on
set(h1,'color',rgb('blueviolet'),'linewidth',2)
set(ax(1),'XColor' , [.3 .3 .3], 'YColor', rgb('blueviolet'), 'LineWidth', 1);
set(ax(1),'YAxisLocation','left','TickDir','out','YMinorTick','off')
set(get(ax(1),'Ylabel'),'String','temperature (°C)')
set(ax(1),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
set(ax(1),'Xdir','reverse','Xlim',[-100 12000])

%set stuff for line plots
set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
set(ax(2),'Xdir','reverse','Xlim',[-100 12000])
set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(get(ax(2),'Ylabel'),'String','z-score');

set(h2,'color',rgb('sandybrown'),'linewidth',2)
set(ax(2),'XColor' , [.3 .3 .3], 'YColor', rgb('sandybrown'), 'LineWidth', 1);
set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(ax(1),'XMinorTick','on','YMinorTick','off', 'YGrid','on')

set(ax(1),'yLim',[-1 0.5])
set(ax(1),'Ytick',-1:0.1:0.5)


%plotting confidence intervals
area_fill(bands(1).time',[calGlobal+allStdGlobal]',[calGlobal-allStdGlobal]',rgb('blueviolet'),rgb('blueviolet'),1,0.2);

area_fill(bands(1).time',[calAnnGlobal+calAnnStdGlobal]',[calAnnGlobal-calAnnStdGlobal]',rgb('black'),rgb('black'),1,0.2);

filled=[[uncalGlobal+uncalStdGlobal]',fliplr([uncalGlobal-uncalStdGlobal]')];
xpoints=[bands(1).time',fliplr(bands(1).time')];

fillhandle=fill(xpoints,filled,rgb('sandybrown'),'Parent',ax(2));%plot the data

set(fillhandle,'EdgeColor',rgb('sandybrown'),'FaceAlpha',0.2,'EdgeAlpha',0.2);%set edge color


% sampledepth subplot
subplot(2,1,2);

%records and sambple bars
[ax,h1,h2] = plotyy(bands(1).time,allSampleGlobal,bands(1).plotEdges,reshape([calSampleGlobal calSampleGlobal]',[],1),@bar,@line); hold on
hold(ax(1))
hold(ax(2))
%other records
plot(ax(2),bands(1).plotEdges,reshape([uncalSampleGlobal uncalSampleGlobal]',[],1), 'color', rgb('sandybrown'), 'linewidth',2);
hold on
plot(ax(2),bands(1).plotEdges,reshape([calAnnSampleGlobal calAnnSampleGlobal]',[],1), 'color', rgb('black'), 'linewidth',2);
hold on

%set stuff for bar graph
set(h1,'facecolor',rgb('Gainsboro'),'edgecolor','none')
set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','left','TickDir','out','YMinorTick','off')
set(ax(1),'yLim',[0 800])
set(ax(2),'yLim',[0 800])
set(ax(2),'Ytick',200:200:800)
set(ax(1),'Ytick',200:200:800)
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
set(h2,'color',rgb('blueviolet'),'linewidth',2)
set(ax(2),'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')