%% Building Figure 6, Temperature 12k, (terrestrial and marine by latitude)
%Cody Routson, July 2019
%% load the things

clear, close all

cols = brewermap(12,'Paired');
load TS.mat
load grid.mat
%% Screen the database

ma = find(strcmp('Temp12k',{TS.paleoData_inCompilation}')); %index the compilation

maa = find(strncmpi('Annual',{TS.interpretation1_seasonalityGeneral}',1)); %index all annual records

mc = find(strncmpi('SummerOnly',{TS.interpretation1_seasonalityGeneral}',1)); %index all SummerOnly records
mi = find(strncmpi('Summer+',{TS.interpretation1_seasonalityGeneral}',1)); %index all SummerOnly records

ms = [mc; mi]; %summer

me = find(strncmpi('WinterOnly',{TS.interpretation1_seasonalityGeneral}',1)); %index all WinterOnly records
mh = find(strncmpi('Winter+',{TS.interpretation1_seasonalityGeneral}',1)); %index all WinterOnly records

mw = [me; mh]; %winter

winter = intersect(mw, ma);
summer = intersect(ms, ma);
annual = intersect(maa, ma);


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
    wint = intersect(winter,banding);
    sum = intersect(summer,banding);
    ann = intersect(annual,banding);
    
    %=======D Bin the records into matrices for summer
    
    for j=1:length(sum) %step though each record in lat band
        %saving the record names in the latitude bands
        Matrix(i).namesSummer{j,1} = TS(sum(j)).dataSetName;
        
        
        sumLatData(j,1) = cell2mat({TS(sum(j)).geo_meanLat});
        sumLonData(j,1) = cell2mat({TS(sum(j)).geo_meanLon});
       
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(sum(j)).age,TS(sum(j)).paleoData_values,binVec);
        if strcmp(TS(sum(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        sumMat(:,j) = bin_mean;
        
    end
    
    clear j
    clear bin_mean
    
    %=======D Bin the records into matrices for the winter
    
    if ~isempty(wint)
        for j=1:length(wint) %step though each record in lat band
            %saving the record names in the latitude bands
            Matrix(i).namesWinter{j,1} = TS(wint(j)).dataSetName;
            
            wintLatData(j,1) = cell2mat({TS(wint(j)).geo_meanLat});
            wintLonData(j,1) = cell2mat({TS(wint(j)).geo_meanLon});
            
            [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(wint(j)).age,TS(wint(j)).paleoData_values,binVec);
            if strcmp(TS(wint(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
                bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
            end
            
            wintMat(:,j) = bin_mean;
            
        end
    else
        wintMat(1:length(BinTime),1) = NaN;
    end
    
    clear j
    clear bin_mean
    
    %=======D Bin the records into matrices for the annual
    
    for j=1:length(ann) %step though each record in lat band
        %saving the record names in the latitude bands
        Matrix(i).namesAnnual{j,1} = TS(ann(j)).dataSetName;
        
        annLatData(j,1) = cell2mat({TS(ann(j)).geo_meanLat});
        annLonData(j,1) = cell2mat({TS(ann(j)).geo_meanLon});
        
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(ann(j)).age,TS(ann(j)).paleoData_values,binVec);
        if strcmp(TS(ann(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        annMat(:,j) = bin_mean;
        
    end
    
    %=======D getting sample depths 
    
    dw = ~isnan(wintMat);
    Matrix(i).wintSampleDepth = nansum(dw,2);
    
    dd = ~isnan(sumMat);
    Matrix(i).sumSampleDepth = nansum(dd,2);
    
    dy = ~isnan(annMat);
    Matrix(i).annSampleDepth = nansum(dy,2);
    
    %=======D saving st*ff 
    
    Matrix(i).wintRecords = (wintMat);
    Matrix(i).sumRecords = (sumMat);
    Matrix(i).annRecords = (annMat);
    Matrix(i).time = binMid;
    Matrix(i).plotEdges = plotEdges;
    
    %=======D Gridding and compositing workflow. 
    
    %SUMMER. Doing if if missing records for a band. 
    if size(sumMat,2)>1 
        Tsunom=zscor_xnan(sumMat); 
        [sumComp,unGridGroups,unGridGroupsLat,unGridGroupsLon,sumGridMean] = gridMat(Tsunom,sumLatData,sumLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
        Matrix(i).sumMedian = sumComp;
    else
        Tsunom=zscor_xnan(dryMat);
        Matrix(i).sumMedian = Tsunom;
    end
    
    
    %WINTER. Doing if if missing records for a band. 
    if size(wintMat,2)>1 
        Twinom=zscor_xnan(wintMat); 
        [winComp,unGridGroups,unGridGroupsLat,unGridGroupsLon,wintGridMean] = gridMat(Twinom,wintLatData,wintLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
        Matrix(i).wintMedian = winComp;
    else
        Twinom=zscor_xnan(wintMat);
        Matrix(i).wintMedian = Twinom;
    end
    
    
     %ANNUAL. Doing "if" if missing records for a band. 
    if size(annMat,2)>1 
        Tannom=zscor_xnan(annMat); 
        [annComp,unGridGroups,unGridGroupsLat,unGridGroupsLon,annGridMean] = gridMat(Tannom,annLatData,annLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
        Matrix(i).annMedian = annComp;
    else
        Tannom=zscor_xnan(annMat);
        Matrix(i).annMedian = Tannom;
    end
    
    %=======D Bootstrap ensembles

    p = bootstrp(500,@nanmedian,sumGridMean');
    Matrix(i).sumEnsemble = p';
    
    q = bootstrp(500,@nanmedian,wintGridMean');
    Matrix(i).wintEnsemble = q';
    
    as = bootstrp(500,@nanmedian,annGridMean');
    Matrix(i).annEnsemble = as';
    

    %clear stuff for next loop
    clear bin_mean
    clear wintMat
    clear sumMat
    clear annMat
    clear band
    clear banding
    clear wintLonData
    clear wintLatData
    clear sumLonData
    clear sumLatData
    clear annLonData
    clear annLatData
    i
end


%% plotting
bands = Matrix;


% Fig5_plot

close all

cc = 1.1;

style_l = {'FontName','Heveltica','FontSize',12,'FontWeight','bold'};
style_i = {'FontName','Heveltica','FontSize',10};


figure

%
FontSize = 16;
FontName = 'Heveltica';

for i = 1:6;

    allSampleDepth = nansum([bands(i).wintSampleDepth,bands(i).sumSampleDepth,bands(i).annSampleDepth],2);
    %0cean. Just use this to scale an axis. Find a better way. 
    ensembleMedian = bands(i).annMedian;
    plotScale1=reshape([ensembleMedian ensembleMedian]',[],1);
    
    % plot sample depths in right panels
    clear sz
    
    %set subplot position and size
    sz = subplot(round(length(bands)),2,i+i);
    sz.Position = sz.Position + [0.05 0 -.05 0]; %set positon and width
   
    %records and sambple bars
    [ax,h1,h2] = plotyy(bands(i).time,allSampleDepth,bands(i).plotEdges,reshape([bands(i).wintSampleDepth bands(i).wintSampleDepth]',[],1),@bar,@line); hold on
    hold(ax(1))
    hold(ax(2))
    %sum records
    plot(ax(2),bands(i).plotEdges,reshape([bands(i).sumSampleDepth bands(i).sumSampleDepth]',[],1), 'color', rgb('firebrick'), 'linewidth',2);
    hold on
    plot(ax(2),bands(i).plotEdges,reshape([bands(i).annSampleDepth bands(i).annSampleDepth]',[],1), 'color', rgb('black'), 'linewidth',2);
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
    
    if max(bands(i).annSampleDepth)>max(bands(i).sumSampleDepth)
        set(ax(2),'yLim',[0 max(bands(i).annSampleDepth)+5])
    else
        set(ax(2),'yLim',[0 max(bands(i).sumSampleDepth)+5])
    end
    
    set(h2,'color',rgb('steelblue'),'linewidth',2)
    set(ax(2),'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
    set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    
    % ticking y at different intervals of sample depth
    if max(bands(i).sumSampleDepth)>max(bands(i).annSampleDepth)
        
        if max(bands(i).sumSampleDepth)>=100
            
            set(ax(2),'Ytick',0:100:roundn(max(bands(i).sumSampleDepth),2))
        else
            set(ax(2),'Ytick',0:5:roundn(max(bands(i).sumSampleDepth),1))
        end
    else
        if max(bands(i).annSampleDepth)>=100
            set(ax(2),'Ytick',0:100:roundn(max(bands(i).sumSampleDepth),2))
        elseif max(bands(i).annSampleDepth)>=50
            set(ax(2),'Ytick',0:20:roundn(max(bands(i).annSampleDepth),1))
        else
            set(ax(2),'Ytick',0:10:roundn(max(bands(i).annSampleDepth),1))
        end
    end

    
    set(get(ax(1),'ylabel'),'String',[bands(i).title]);
    set(get(ax(1),'Ylabel'),'FontName','Heveltica','FontSize',12,'color',[0 0 0]),
    
    
    % plot records in left panels
    clear sx
    %set subplot position and size
    sx = subplot(round(length(bands)),2,i+i-1);
    sx.Position = sx.Position + [0.0 0 .05 0];
    
    %sum composites
    if max(bands(i).sumSampleDepth)>=5
        ha = plot(bands(i).time,bands(i).sumMedian,'linewidth',2,'color',rgb('firebrick'));
        area_fill(bands(i).time',[bands(i).sumMedian+nanstd(bands(i).sumEnsemble,0,2)]',[bands(i).sumMedian-nanstd(bands(i).sumEnsemble,0,2)]',rgb('firebrick'),rgb('firebrick'),1,0.2);
    end
    hold on
    %one composites
    if max(bands(i).wintSampleDepth)>=5
        hb = plot(bands(i).time,bands(i).wintMedian,'linewidth',2,'color',rgb('steelblue'));
        area_fill(bands(i).time',[bands(i).wintMedian+nanstd(bands(i).wintEnsemble,0,2)]',[bands(i).wintMedian-nanstd(bands(i).wintEnsemble,0,2)]',rgb('steelblue'),rgb('steelblue'),1,0.2);
    end
    hold on
    %one composites
    if max(bands(i).annSampleDepth)>=5
        hc = plot(bands(i).time,bands(i).annMedian,'linewidth',2,'color',rgb('black'));
        area_fill(bands(i).time',[bands(i).annMedian+nanstd(bands(i).annEnsemble,0,2)]',[bands(i).annMedian-nanstd(bands(i).annEnsemble,0,2)]',rgb('gray'),rgb('gray'),1,0.2);
    end
    %making things look a little better
    hold on
    set(gca,'Xdir','reverse','Xlim',[-50 12000])
    set(get(gca,'Ylabel'),'String','z-score')
    set(gca,'FontName',FontName,'FontSize', round(FontSize*0.71));
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
    set(gca,'yLim',[nanmean(plotScale1)-1.5 nanmean(plotScale1)+1.5])

    
end

