%% Building Figure 5, Temperature 12k. By season and latitude
%Cody Routson, July 2019
%% load the things

clear
cols = brewermap(12,'Paired');
load('Temp12k_v1_0_0.mat','TS')
load grid.mat
%% Screen the database
ma = find(strcmp('Temp12k',{TS.paleoData_inCompilation}')); %index the compilation

mn = find(strncmpi('MarineSediment',{TS.archiveType}',1)); %index all marine archives

mb = find(strncmpi('Annual',{TS.interpretation1_seasonalityGeneral}',1)); %index all annual records
mc = find(strncmpi('SummerOnly',{TS.interpretation1_seasonalityGeneral}',1)); %index all SummerOnly records
me = find(strncmpi('WinterOnly',{TS.interpretation1_seasonalityGeneral}',1)); %index all WinterOnly records

mf = [mb; mc; me]; %index of annual, summer and winter

mg = intersect(mf,ma); %intersect seasons with compilation for a USE IN index

marine = intersect(mn,mg); % intersect marine with the USE IN index to find marine index
terrest = setdiff(mg,marine); %set difference between USE IN and marine to find terrestrial index

%% grabbing all latitudes 

latitude={TS.geo_latitude};
lat = (cell2mat(latitude))';

%% preallocating the latitudinal bands

Matrix(6).index = find(lat>=-90 & lat<-60);
Matrix(6).title = '90°-60°S';

Matrix(5).index = find(lat>=-60 & lat<-30);
Matrix(5).title = '60°-30°S';

Matrix(4).index = find(lat>=-30 & lat<0);
Matrix(4).title = '30°-0°S';

Matrix(3).index = find(lat>=0 & lat<30);
Matrix(3).title = '0°-30°N';

Matrix(2).index = find(lat>=30 & lat<60);
Matrix(2).title = '30°-60°N';

Matrix(1).index = find(lat>=60 & lat<90);
Matrix(1).title = '60°-90°N';

%% latitude steps for bands
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
    dirt = intersect(terrest,banding);
    wet = intersect(marine,banding);
    
    %%=======D Bin the records into matrices for marine
    
    for j=1:length(wet) %step though each ocean record in lat band
        %saving the record names in the latitude bands
        
        
        wetLatData(j,1) = cell2mat({TS(wet(j)).geo_meanLat});
        wetLonData(j,1) = cell2mat({TS(wet(j)).geo_meanLon});
        
        Matrix(i).namesMarine{j,1} = TS(wet(j)).dataSetName;
       
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(wet(j)).age,TS(wet(j)).paleoData_values,binVec);
        if strcmp(TS(wet(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        wetMat(:,j) = bin_mean;
        
    end
    
    clear j
    clear bin_mean
    
    
    %=======D Bin the records into matrices for terrestrial
    
    
    for j=1:length(dirt) %step though each terrestrial record in lat band
        %saving the record names in the latitude bands
        Matrix(i).namesTerrest{j,1} = TS(dirt(j)).dataSetName;
        
        dryLatData(j,1) = cell2mat({TS(dirt(j)).geo_meanLat});
        dryLonData(j,1) = cell2mat({TS(dirt(j)).geo_meanLon});
        
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(dirt(j)).age,TS(dirt(j)).paleoData_values,binVec);
        if strcmp(TS(dirt(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        dryMat(:,j) = bin_mean;
        
    end
    
    %=======D getting sample depths for marine and terrestrial
    
    dw = ~isnan(wetMat);
    Matrix(i).oceanSampleDepth = sum(dw,2);
    
    dd = ~isnan(dryMat);
    Matrix(i).terrestSampleDepth = sum(dd,2);
    
    %=======D saving st*ff 
    
    Matrix(i).oceanRecords = (wetMat);
    Matrix(i).terrestRecords = (dryMat);
    Matrix(i).time = binMid;
    Matrix(i).plotEdges = plotEdges;

    %=======D Gridding and compositing workflow. 
    
    %uncalibrated records
    oMean=nanmean(wetMat(binMid>normStart & binMid<normEnd,:),1); %using the mean for some period to compute anomalies %Default is to normalize over the whole 12000 year interval. 
    oMat=repmat(oMean,length(BinTime),1); % turn it into a matrix for subtraction
    Tanom=zscor_xnan(wetMat);  % subract out the average
    [wetComp,gridGroups,gridGroupsLat,gridGroupsLon,gridMean] = gridMat(Tanom,wetLatData,wetLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
    Matrix(i).oceanMedian = wetComp;
    
    %uncalibrated records. Doing if if missing records for a band. Not an issue here
    if size(dryMat,2)>1 %~isempty(uncal)
        Tunom=zscor_xnan(dryMat); % zscores for uncalibrated rcords
        [uncalComp,unGridGroups,unGridGroupsLat,unGridGroupsLon,unGridMean] = gridMat(Tunom,dryLatData,dryLonData,latgrid,longrid); %GRIDD run gridMat and collect total mean into a matrix
        Matrix(i).terrestMedian = uncalComp;
    else
        Tunom=zscor_xnan(dryMat);
        Matrix(i).terrestMedian = Tunom;
    end
     
    %=======D Bootstrap ensembles

    m = bootstrp(500,@nanmedian,gridMean');
    
    Matrix(i).oceanEnsemble = m';
    
    p = bootstrp(500,@nanmedian,unGridMean');
    Matrix(i).terrestEnsemble = p';
    
    
    %clear stuff for next loop
    clear bin_mean
    clear wetMat
    clear dryMat
    clear band
    clear banding
    clear wetLatData
    clear wetLonData
    clear dryLatData
    clear dryLonData
    i
end


%% scheming and plotting
bands = Matrix;

cc = 1.1;

style_l = {'FontName','Heveltica','FontSize',12,'FontWeight','bold'};
style_i = {'FontName','Heveltica','FontSize',10};

fig('Fig5')


%
FontSize = 16;
FontName = 'Heveltica';

for i = 1:6;
    
    
    allSampleDepth = sum([bands(i).oceanSampleDepth,bands(i).terrestSampleDepth],2);
    %0cean. Just using this to scale an axis. Find a better way. 
    ensembleMedian = bands(i).oceanMedian;
    plotScale1=reshape([ensembleMedian ensembleMedian]',[],1);
    
    %Terestrial. Don't end up using this
    ensembleMedian2 = bands(i).terrestMedian;
    plotScale2=reshape([ensembleMedian2 ensembleMedian2]',[],1);
    
    
    % plot sample depths in right panels
    clear sz
    
    %set subplot position and size
    sz = subplot(round(length(bands)),2,i+i);
    sz.Position = sz.Position + [0.05 0 -.05 0]; %set positon and width

    %marine records and sambple bars
    [ax,h1,h2] = plotyy(bands(i).time,allSampleDepth,bands(i).plotEdges,reshape([bands(i).oceanSampleDepth bands(i).oceanSampleDepth]',[],1),@bar,@line); hold on %% 'color', rgb('darkcyan'), 'linewidth',2);
    hold(ax(1))
    hold(ax(2))
    
    %terrestrial records
    plot(ax(2),bands(i).plotEdges,reshape([bands(i).terrestSampleDepth bands(i).terrestSampleDepth]',[],1), 'color', rgb('peru'), 'linewidth',2);
    
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
    
    if max(bands(i).terrestSampleDepth)>max(bands(i).oceanSampleDepth)
        set(ax(2),'yLim',[0 max(bands(i).terrestSampleDepth)+5])
    else
        set(ax(2),'yLim',[0 max(bands(i).oceanSampleDepth)+5])
    end
    
    set(h2,'color',rgb('darkcyan'),'linewidth',2)
    set(ax(2),'XColor' , [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
    set(ax(2),'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    
    
    if max(bands(i).terrestSampleDepth)>max(bands(i).oceanSampleDepth)
        
        if max(bands(i).terrestSampleDepth)>=100 && max(bands(i).terrestSampleDepth)<400
            
            set(ax(2),'Ytick',0:100:roundn(max(bands(i).terrestSampleDepth),2))
        elseif max(bands(i).terrestSampleDepth)>=400
            set(ax(2),'Ytick',0:200:roundn(max(bands(i).terrestSampleDepth),2))
        else
            set(ax(2),'Ytick',0:5:roundn(max(bands(i).terrestSampleDepth),1))
        end
    else
        if max(bands(i).oceanSampleDepth)>=50 
            set(ax(2),'Ytick',0:20:roundn(max(bands(i).oceanSampleDepth),1))
        else
            set(ax(2),'Ytick',0:10:roundn(max(bands(i).oceanSampleDepth),1))
        end
    end

%     
    set(get(ax(1),'ylabel'),'String',[bands(i).title]);
    set(get(ax(1),'ylabel'),'FontName','Heveltica','FontSize',12,'color',[0 0 0]),
     
    % plot records in left panels
    clear sx
    %set subplot position and size
    sx = subplot(round(length(bands)),2,i+i-1);
    sx.Position = sx.Position + [0.0 0 .05 0];
    
    %terrestrial composites
    if max(bands(i).terrestSampleDepth)>=5
        ha = plot(bands(i).time,bands(i).terrestMedian,'linewidth',2,'color',rgb('peru'));
        area_fill(bands(i).time',[bands(i).terrestMedian+nanstd(bands(i).terrestEnsemble,0,2)]',[bands(i).terrestMedian-nanstd(bands(i).terrestEnsemble,0,2)]',rgb('peru'),rgb('peru'),1,0.2);
        hold on
    end
    %marine composites
    if max(bands(i).oceanSampleDepth)>=5
        hb = plot(bands(i).time,bands(i).oceanMedian,'linewidth',2,'color',rgb('darkcyan'));
        area_fill(bands(i).time',[bands(i).oceanMedian+nanstd(bands(i).oceanEnsemble,0,2)]',[bands(i).oceanMedian-nanstd(bands(i).oceanEnsemble,0,2)]',rgb('darkcyan'),rgb('darkcyan'),1,0.2);
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

%% getting the data to save for later
for j = 1:6
    BA.a(:,j) = bands(j).terrestMedian;
    BA.b(:,j) = bands(j).oceanMedian;
    
    BA.asd(:,j) = nanstd(bands(j).terrestEnsemble,0,2);
    BA.bsd(:,j) = nanstd(bands(j).oceanEnsemble,0,2);
    
    BA.aa(:,j) = bands(j).terrestSampleDepth;
    BA.bb(:,j) = bands(j).oceanSampleDepth;
end
