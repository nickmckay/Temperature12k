%% Building Figure 7, Temperature 12k. High and low resolution by latitude. 
%Cody Routson, July 2019
%% load the things

clear, close all

cols = brewermap(12,'Paired');
load TS.mat
load grid.mat

%% Screen the database
maxRes = 400; %SET RESOLUTION (MAXIMUM INTERVAL BETWEEN SAMPLES)
minRes = 100;


allYears = {TS.age};

age10k = cellfun(@(x) x(x<12000), allYears,'UniformOutput',0);
emp =find(cellfun(@isempty,age10k));
age10k(emp)={NaN};

medRes = cellfun(@(x) median(abs(diff(x))),age10k);


%% Screening for high and low resolution records in Temp12k
mch = find(strcmp('Temp12k',{TS.paleoData_inCompilation}') & (medRes <= minRes)');

mcf = find(strcmp('Temp12k',{TS.paleoData_inCompilation}') & (medRes > minRes)');

mb = find(strncmpi('annual',{TS.interpretation1_seasonalityGeneral}',1)); %index all annual records
mc = find(strncmpi('summerOnly',{TS.interpretation1_seasonalityGeneral}',1)); %index all SummerOnly records
me = find(strncmpi('winterOnly',{TS.interpretation1_seasonalityGeneral}',1)); %index all WinterOnly records

mf = [mb; mc; me]; %index of annual, summer and winter

high = intersect(mf,mch);
low = intersect(mf, mcf);

%% grabbing all latitudes 

latitude={TS.geo_latitude};
lat = (cell2mat(latitude))';

%% preallocating the latitudinal bands

Matrix(6).index = find(lat>=-90 & lat<-60);
Matrix(6).title = '90°-60°S';

Matrix(5).index = find(lat>=-60 & lat<-30);
Matrix(5).title = '60°-30°S';

Matrix(4).index = find(lat>=-30 & lat<05);
Matrix(4).title = '30°-0°S';

Matrix(3).index = find(lat>=0 & lat<30);
Matrix(3).title = '0°-30°N';

Matrix(2).index = find(lat>=30 & lat<60);
Matrix(2).title = '30°-60°N';

Matrix(1).index = find(lat>=60 & lat<90);
Matrix(1).title = '60°-90°N';



%% latitude steps
latstep = 90:-30:-90;


%% preallocate Bins
 %binning at lower resolution to calculate sample depts
binStep=500; %SET BINSTEP
binEdges=[0 12000]; %SET PERIOD OF ANALYSIS (yr BP)
binVec=min(binEdges):binStep:max(binEdges);
binMid=binVec(1:end-1)+binStep/2;
binVec = binVec';
binMid = binMid';
plotEdges=reshape([binVec(1:end-1) binVec(2:end)]',[],1);

%% binning for high resolution
binStep2=100; %SET BINSTEP
binEdges2=[0 12000]; %SET PERIOD OF ANALYSIS (yr BP)
binVec2=min(binEdges2):binStep2:max(binEdges2);
binMid2=binVec2(1:end-1)+binStep2/2;
binVec2 = binVec2';
binMid2 = binMid2';
plotEdges2=reshape([binVec2(1:end-1) binVec2(2:end)]',[],1);

%% set the equal areas for griddng 
normStart = 0; %SET START OF NORMALIZATION PERIOD (yr BP)
normEnd = 12000; %SET END OF NORMALIZATION PERIOD (yr BP)

latgrid = grid.latgrid;
longrid = grid.longrid;
% [latgrid, longrid, regbound] = equalGridDeg(4000, -90);

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
    lowRes = intersect(low,banding);
    highRes = intersect(high,banding);
    
 %=======D Bin the records into matrices for the HIGHRES
 
    for j=1:length(highRes) %step though each record in lat band
        %saving the record names in the latitude bands
        Matrix(i).namesHighRes{j,1} = TS(highRes(j)).dataSetName;
        
        highLatData(j,1) = cell2mat({TS(highRes(j)).geo_meanLat});
        highLonData(j,1) = cell2mat({TS(highRes(j)).geo_meanLon});
       
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(highRes(j)).age,TS(highRes(j)).paleoData_values,binVec2);
        [bin_mean2, BinTime2, sem2, bin_sum2] =  bin_x(TS(highRes(j)).age,TS(highRes(j)).paleoData_values,binVec); %binning at lower resolution to calculate total sample depths
        
        if strcmp(TS(highRes(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        highMat(:,j) = (bin_mean);
        highMat2(:,j) = (bin_mean2);
        
    end
    
    clear j
    clear bin_mean
    clear bin_mean2
    
    %=======D Bin the records into matrices for the LOWRES
    
    if ~isempty(lowRes)
    for j=1:length(lowRes) %step though each record in lat band
        %saving the record names in the latitude bands
        Matrix(i).namesLowRes{j,1} = TS(lowRes(j)).dataSetName;
        
        lowLatData(j,1) = cell2mat({TS(lowRes(j)).geo_meanLat});
        lowLonData(j,1) = cell2mat({TS(lowRes(j)).geo_meanLon});
        
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(lowRes(j)).age,TS(lowRes(j)).paleoData_values,binVec);
        if strcmp(TS(lowRes(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        lowMat(:,j) = nanzscore(bin_mean);
        
    end
    else
        lowMat(1:24,1) = NaN;
    end
    
    
    %=======D getting sample depths 
    
    dw = ~isnan(highMat);
    Matrix(i).highSampleDepth = sum(dw,2);
    
    dh = ~isnan(highMat2); %Sample depth of high resolution records at low resolution binning
    Matrix(i).highSDP = sum(dh,2);
    
    dd = ~isnan(lowMat);
    Matrix(i).lowSampleDepth = sum(dd,2);
    
    %=======D saving st*ff 
    
    Matrix(i).highRecords = (highMat);
    Matrix(i).lowRecords = (lowMat);
    Matrix(i).time = binMid;
    Matrix(i).timeHigh = binMid2;
    Matrix(i).plotEdgesHigh = plotEdges2;
    Matrix(i).plotEdges = plotEdges;
    

    %=======D Gridding and compositing workflow. 
    
    
      %HIGH RESOLUTION. Doing "if" if missing records for a band. 
    if size(highMat,2)>1 
        Thinom=zscor_xnan(highMat); 
        [highComp,unGridGroups,unGridGroupsLat,unGridGroupsLon,highGridMean] = gridMat(Thinom,highLatData,highLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
        Matrix(i).highMedian = highComp;
    else
        Thinom=zscor_xnan(highMat);
        Matrix(i).highMedian = Thinom;
    end
    
    
    %LOW RESOLUTION. Doing "if" if missing records for a band.
    if size(lowMat,2)>1 
        Tlonom=zscor_xnan(lowMat); 
        [lowComp,unGridGroups,unGridGroupsLat,unGridGroupsLon,lowGridMean] = gridMat(Tlonom,lowLatData,lowLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
        Matrix(i).lowMedian = lowComp;
    else
        Tlonom=zscor_xnan(lowMat);
        Matrix(i).lowMedian = Tlonom;
    end
    
    
    %=======D Bootstrap ensembles

    p = bootstrp(500,@nanmedian,highGridMean');
    Matrix(i).highEnsemble = p';
    
    q = bootstrp(500,@nanmedian,lowGridMean');
    Matrix(i).lowEnsemble = q';
    
    
    %clear stuff for next loop
    clear bin_mean
    clear highMat
    clear highMat2
    clear lowMat
    clear band
    clear banding
    clear latData
    clear lonData
    clear dd
    clear dh
    clear lowLonData
    clear lowLatData
    clear highLonData
    clear highLatData
    i
end


%% Sending latitudinal bands away for plotting
bands = Matrix;


% Fig6_plot
cc = 1.1;



cols = [cols; 0.2157    0.4941    0.7216];

style_l = {'FontName','Heveltica','FontSize',12,'FontWeight','bold'};
style_i = {'FontName','Heveltica','FontSize',10};


figure

%
FontSize = 16;
FontName = 'Heveltica';

for i = 1:6;

    allSampleDepth = nansum([bands(i).highSDP,bands(i).lowSampleDepth],2); 
    %Just use this to scale an axis. Find a better way. 
    ensembleMedian = bands(i).lowMedian;
    plotScale1=reshape([ensembleMedian ensembleMedian]',[],1);
    
    % plot sample depths in right panels
    clear sz
    
    %set subplot position and size
    sz = subplot(round(length(bands)),2,i+i);
    sz.Position = sz.Position + [0.05 0 -.05 0]; %set positon and width
   

    %records and sambple bars
    [ax,h1,h2] = plotyy(bands(i).time,allSampleDepth,bands(i).plotEdgesHigh,reshape([bands(i).highSampleDepth bands(i).highSampleDepth]',[],1),@bar,@line); hold on 
    hold(ax(1))
    hold(ax(2))
    %other records
    plot(ax(2),bands(i).plotEdges,reshape([bands(i).lowSampleDepth bands(i).lowSampleDepth]',[],1), 'color', rgb('midnightblue'), 'linewidth',2);
    hold on
    
    %set stuff for bar graph
    set(h1,'facecolor',rgb('Gainsboro'),'edgecolor','none')
    set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','left','TickDir','out','YMinorTick','off')
    set(ax(1),'yLim',[0 800])
    set(ax(1),'Ytick',200:200:800)
    set(ax(1),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
    set(ax(1),'Xdir','reverse','Xlim',[-100 12000])
    
    %set stuff for line plots
    set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
    set(ax(2),'Xdir','reverse','Xlim',[-100 12000])
    set(get(ax(2),'Ylabel'),'String','#records');
    set(get(ax(2),'Ylabel'),style_l{:}),
    set(gcf,'CurrentAxes',ax(2))
    
    if max(bands(i).lowSampleDepth)>max(bands(i).highSampleDepth)
        set(ax(2),'yLim',[0 max(bands(i).lowSampleDepth)+5])
    else
        set(ax(2),'yLim',[0 max(bands(i).highSampleDepth)+5])
    end
    
    set(h2,'color',rgb('dimgray'),'linewidth',2)
    set(ax(2),'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
    set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    
    if max(bands(i).lowSampleDepth)>max(bands(i).highSampleDepth)
        
        if max(bands(i).lowSampleDepth)>=100
            
            set(ax(2),'Ytick',0:100:roundn(max(bands(i).lowSampleDepth),2))
        elseif max(bands(i).lowSampleDepth)>=50
            set(ax(2),'Ytick',0:20:roundn(max(bands(i).lowSampleDepth),1))
            
        else
            set(ax(2),'Ytick',0:10:roundn(max(bands(i).lowSampleDepth),1))
        end
    else
        set(ax(2),'Ytick',0:5:round(max(bands(i).highSampleDepth)/5)*5)
    end

    set(get(ax(1),'ylabel'),'String',[bands(i).title]);
    set(get(ax(1),'Ylabel'),'FontName','Heveltica','FontSize',12,'color',[0 0 0]),
    
    
    % plot records in left panels
    clear sx
    %set subplot position and size
    sx = subplot(round(length(bands)),2,i+i-1);
    sx.Position = sx.Position + [0.0 0 .05 0];
    
    %other composites
    if max(bands(i).lowSampleDepth)>=5
        ha = plot(bands(i).time,bands(i).lowMedian,'linewidth',2,'color',rgb('midnightblue'));
        area_fill(bands(i).time',[bands(i).lowMedian+nanstd(bands(i).lowEnsemble,0,2)]',[bands(i).lowMedian-nanstd(bands(i).lowEnsemble,0,2)]',rgb('midnightblue'),rgb('midnightblue'),1,0.2);
    end
    hold on
    %one composites
    if max(bands(i).highSampleDepth)>=5
        hb = plot(bands(i).timeHigh,bands(i).highMedian,'linewidth',2,'color',rgb('dimgray'));
        area_fill(bands(i).timeHigh',[bands(i).highMedian+nanstd(bands(i).highEnsemble,0,2)]',[bands(i).highMedian-nanstd(bands(i).highEnsemble,0,2)]',rgb('dimgray'),rgb('dimgray'),1,0.2);
    end
    %making things look a little better
    hold on
    set(gca,'Xdir','reverse','Xlim',[-50 12000])
    set(get(gca,'Ylabel'),'String','z-score')
    set(gca,'FontName',FontName,'FontSize', round(FontSize*0.71));
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
    set(gca,'yLim',[nanmean(bands(i).highMedian)-1.5 nanmean(bands(i).highMedian)+1.5])

    clear allSampleDepth
end




