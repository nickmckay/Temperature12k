%% Building global composite using calibrated records and 3000-5000 year interval of overlap, Temperature 12k. 
%Cody Routson, November 2019
%% load the things

clear, close all

cols = brewermap(12,'Paired');
load TS.mat
load PAGES_multiMethodMeadian.txt;
load grid.mat;

%% Interval over which to subtract the mean (calibrated records) and normalize (uncalibrated records)
normStart = 3000; %SET START OF NORMALIZATION PERIOD (yr BP)
normEnd = 5000; %SET END OF NORMALIZATION PERIOD (yr BP)
binStep=100; %SET BINSTEP
binEdges=[-50 11950]; %SET PERIOD OF ANALYSIS (yr BP)
loops = 500;

%% Screen the database

mn = find(strncmpi('degC',{TS.paleoData_units}',1)); %index all calibrated
ma = find(strcmp('Temp12k',{TS.paleoData_inCompilation}')); %index Temp12k
mb = find(strncmpi('annual',{TS.interpretation1_seasonalityGeneral}',7)); %index all annual records
mc = find(strncmp('summerOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all SummerOnly records
me = find(strncmp('winterOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all WinterOnly records

mf = [mb; mc; me]; %index of annual, summer and winter

mg = intersect(mf,ma); %seasons and Temp12k index

gm = intersect(mb,ma); %annual and Temp12k index

zog = intersect(mn,mg); % finding calibrated index


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

binVec=min(binEdges):binStep:max(binEdges);
binMid=binVec(1:end-1)+binStep/2;
binVec = binVec';
binMid = binMid';

plotEdges=reshape([binVec(1:end-1) binVec(2:end)]',[],1);

%% equal area latitude and longitude grid
latgrid = grid.latgrid;
longrid = grid.longrid;

%% calculate weight
w(1,1) = (sind(90)-sind(60))/2;
w(1,2) = (sind(60)-sind(30))/2;
w(1,3) = (sind(30)-sind(0))/2;
w(1,4) = (sind(30)-sind(0))/2;
w(1,5) = (sind(60)-sind(30))/2;
w(1,6) = (sind(90)-sind(60))/2;

w1 = repmat(w, 24,1);

%% Loop though each latitude 30° band, Index records. Bin records. Add error. Calculate composites. Repeat over bootstrapped iterations (loops).  
for ii = 1:loops;
    for i = 1:6; %step through each 30° latitude band
        
        % Start the playing in a band (index all rows in latitudinal range)
        band = [latstep(i+1), latstep(i)];
        banding = find(cell2mat(latitude)>band(1,1) & cell2mat(latitude)< band(1,2)); %% specify the latitude cuttoff
        banding = banding';
        
        % index of records within latitude band (intersecting the records we
        % want in the latitudinal intervals
        c = intersect(zog,banding); %calibrated
        
        %=======D Bin the records into matrices
        
        for j=1:length(c) %step though each record in lat band
            %saving the record names in the latitude bands
            Matrix(i).namesCal{j,1} = TS(c(j)).dataSetName;
            
            %put site level lat and lon data into a matrix for gridding later
            calLatData(j,1) = cell2mat({TS(c(j)).geo_meanLat});
            calLonData(j,1) = cell2mat({TS(c(j)).geo_meanLon});
            
           er = TS(c(j)).paleoData_temperature12kUncertainty;
           if strncmp(er,'NA',2) == 1
               er = 1.5;
           end
           if isempty(er) == 1
               er = 1.5;
           end
            
            %bin the records
            if ii == 1; %for first iteration use the data without adding any uncertainty
                [bin_mean1, BinTime, sem, bin_sum] =  bin_x(TS(c(j)).age,TS(c(j)).paleoData_values,binVec);
                calMat1(:,j) = bin_mean1;
            else %for other xxx iterations, add age and proxy uncertainties from a normal distrubution
                [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(c(j)).age*normrnd(1,0.05),TS(c(j)).paleoData_values+normrnd(0,er,size(TS(c(j)).paleoData_values,1),1),binVec);
                calMat(:,j) = bin_mean;
            end
            
        end
        
        clear j
        clear bin_mean
        clear bin_mean1
        
        if ii == 1
            %=======D saving st*ff
            
            Matrix(i).calRecords = (calMat1);
            Matrix(i).time = binMid;
            Matrix(i).plotEdges = plotEdges;

            %=======D Gridding and compositing workflow.
            
            % making anomalies gridding calibrated records
            oMean=nanmean(calMat1(binMid>normStart & binMid<normEnd,:),1); %using the mean for some period to subtract from each record %Default is to subtract the mean of the whole 12000 year interval.
            oMat=repmat(oMean,length(BinTime),1); % turn it into a matrix for subtraction
            Tanom=calMat1-oMat; % subract out the average
            [calComp,gridGroups,gridGroupsLat,gridGroupsLon,gridMean] = gridMat(Tanom,calLatData,calLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
            
            Matrix(i).calMedian = calComp;
            Matrix(i).calEnsemble(:,ii) = calComp;
        
            %=======D getting sample depths
            
            dw = ~isnan(Tanom);
            Matrix(i).calSampleDepth = sum(dw,2);
            Matrix(i).sampleN = dw;
            
            %=======D getting number of records for each band
            
            tes = sum(Matrix(i).sampleN,1);
            ab = find(tes>0);
            Matrix(i).recordN = length(ab);

        else
            
            %=======D Gridding and compositing workflow.
            
            % making anomalies gridding calibrated records
            oMean=nanmean(calMat(binMid>normStart & binMid<normEnd,:),1); %using the mean for some period to subtract from each record %Default is to subtract the mean of the whole 12000 year interval.
            oMat=repmat(oMean,length(BinTime),1); % turn it into a matrix for subtraction
            Tanom=calMat-oMat; % subract out the average
            [calComp,gridGroups,gridGroupsLat,gridGroupsLon,gridMean] = gridMat(Tanom,calLatData,calLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
            
            Matrix(i).calEnsemble(:,ii) = calComp;
            
        end
        
        %clear st*ff for next loop
        clear bin_mean
        clear calMat
        clear calMat1
        clear calAnnMat
        clear uncalMat
        clear band
        clear banding
        clear calLatData
        clear calLonData
        clear calComp

    end
    
    % create a global composite for each loop
    %multiply weights
    for l = 1:6
        cy(:,l) = w1(1,l)*Matrix(l).calEnsemble(:,ii);
    end
    
    globalComp.medianEnsemble(:,ii) = sum(cy,2);
    
    ii
end
%%
bands = Matrix;

%% Looping through the latitudinal bands for global sample depths
for k=1:6 
    calSampleDepth(:,k) = bands(k).calSampleDepth;
    bands(k).standardDev = nanstd(bands(k).calEnsemble,0,2);
    recordN(:,k) = bands(k).recordN;
end

calSampleGlobal = nansum(calSampleDepth,2);
calSmoothGlobal = nanmedian(globalComp.medianEnsemble,2);
allStdGlobal = nanstd(globalComp.medianEnsemble,0,2);

%% subtracting the mean from the pages data
p2k.time = PAGES_multiMethodMeadian(:,1);
p2k.tmp = PAGES_multiMethodMeadian(:,2);

gd = find(p2k.time>=500 & p2k.time<=1500);

gt = nanmean(p2k.tmp(gd));
p2k.tmp = p2k.tmp-gt;

%bin the pages
[bin_mean2k, BinTime2k, sem, bin_sum] =  bin_x(p2k.time,p2k.tmp,binVec);

%get 1850-1900 registered to zero
gt = find(p2k.time>=50 & p2k.time<=100);
gx = nanmean(p2k.tmp(gt));
bin_mean2k = bin_mean2k-gx;
p2k.tmp = p2k.tmp-gx;

%% global composite without scale structure to save

globalComp.recordN = nansum(recordN,2);
globalComp.time = bands(1).time;
globalComp.globalComposite = globalComp.medianEnsemble(:,1);
globalComp.sampleDepth = calSampleGlobal;
globalComp.ensembleMedian = calSmoothGlobal;
globalComp.stDev = allStdGlobal;
globalComp.p2kTime = p2k.time;
globalComp.p2kTemp = p2k.tmp;
globalComp.p2k200Bin = bin_mean2k;
globalComp.p2k200Time = BinTime2k;
globalComp.mean1800to1900 = calSmoothGlobal(2);

xp = calSmoothGlobal(2);
%% save some j*nk

save latBands30Deg.mat bands

save globalComp.mat globalComp

%% Plot
style_l = {'FontName','Heveltica','FontSize',12,'FontWeight','bold'};
figure(1)
clf

%timeseries subplot
ci_archs.quant  = quantile(globalComp.medianEnsemble,[0.025 0.975],2);
ciLo = ci_archs.quant(:,1)-xp;% reshape([ci_archs.quant(:,1) ci_archs.quant(:,1)]',[],1);
ciHi = ci_archs.quant(:,2)-xp;%reshape([ci_archs.quant(:,2) ci_archs.quant(:,2)]',[],1);

ylims       = [floor(1.1*min(ciLo)), ceil(1.1*max(ciHi))];

[ax,h1,h2] = plotyy(bands(1).time,calSampleGlobal ,bands(1).time,calSmoothGlobal-xp,@bar,@line); hold on
hold(ax(2))
hold on

hold on

set(h1,'facecolor',rgb('Gainsboro'),'edgecolor','none')
set(ax(2),'XColor' , [.3 .3 .3], 'YColor', rgb('steelblue'), 'LineWidth', 1);
set(ax(2),'YAxisLocation','left','TickDir','out','YMinorTick','off')
set(get(ax(2),'Ylabel'),'String','temperature (°C)')
set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
set(ax(2),'Xdir','reverse','Xlim',[-100 12000])

%set stuff
set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right','TickDir','out','YMinorTick','off')
set(ax(1),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
set(ax(1),'yLim',[0 600])
set(ax(1),'Xdir','reverse','Xlim',[-100 12000])
set(ax(1),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(get(ax(1),'Ylabel'),'String','# records');

set(h2,'color',rgb('steelblue'),'linewidth',2)
set(ax(1),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
set(ax(2),'XMinorTick','on','YMinorTick','off', 'YGrid','on')

set(ax(2),'Ytick',-1:0.1:0.75)

set(gcf,'CurrentAxes',ax(2))

hci = area_fill(bands(1).time',ciLo',ciHi',rgb('steelblue'),rgb('steelblue'),0.1);


%%
figure(2)
clf
hold on
for i=1:6
    txt = bands(i).title;
    plot(bands(i).time,bands(i).calMedian,'DisplayName',txt)
end
hold off
legend show
legend boxoff
set(gca,'Xdir','reverse','Xlim',[-100 12000])
xlabel('year (BP)')
ylabel('temperature (°C)')

