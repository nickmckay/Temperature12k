%% proxy plot, Temperature 12k (Figure 4)
%Cody Routson, Oct 2019

clear

load('Temp12k_v1_0_0.mat','TS')

binStep=500; %SET BINSTEP
binEdges=[0 12000]; %SET PERIOD OF ANALYSIS (yr BP)
normStart = 0; %SET START OF NORMALIZATION PERIOD (yr BP)
normEnd = 12000; %SET END OF NORMALIZATION PERIOD (yr BP)
%% Indexing the data
ma = find(strcmp('Temp12k',{TS.paleoData_inCompilation}')); %index the compilation

mb = find(strncmp('annual',{TS.interpretation1_seasonalityGeneral}',7)); %index all annual records
mc = find(strncmp('summerOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all SummerOnly records
me = find(strncmp('winterOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all WinterOnly records
mf = [mb; mc; me]; %index of annual, summer and winter

mg = intersect(mf,ma); %intersect seasons with compilation for a USE IN index
%% proxy specific indexes
mgca = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'Mg/Ca'),TS));
alk = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'alkenone'),TS));
pollen = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'pollen'),TS));
chironomid = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'chironomid'),TS));
biophysical = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'biophysical'),TS));
biomarker = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'other biomarker'),TS));
isotope = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'isotope'),TS));
microfossil = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'other microfossil'),TS));
otherIce = find(arrayfun(@(x)strcmp(x.paleoData_proxyGeneral,'other ice'),TS));

%% preallocating stuff for the plotting into a structure

Matrix(8).index =  intersect(mg,microfossil);
Matrix(8).title = 'other microfossil';
Matrix(8).color = [128/255 54/255 168/255]; 

Matrix(7).index =  intersect(mg,isotope);
Matrix(7).title = 'isotope';
Matrix(7).color = [236/255 156/255 46/255];  

Matrix(6).index =  intersect(mg,biomarker);
Matrix(6).title = 'other biomarker';
Matrix(6).color = [33/255 52/255 219/255]; 

Matrix(5).index = intersect(mg,biophysical);
Matrix(5).title = 'biophysical';
Matrix(5).color = [53/255 165/255 197/255];

Matrix(4).index = intersect(mg,alk);
Matrix(4).title = 'alkenone';
Matrix(4).color = [221/255 84/255 46/255];

Matrix(3).index = intersect(mg,pollen);
Matrix(3).title = 'pollen';
Matrix(3).color = [50/255 127/255 81/255]; 

Matrix(2).index = intersect(mg,mgca);
Matrix(2).title = 'Mg/Ca';
Matrix(2).color = [8/255 46/255 114/255]; 

Matrix(1).index = intersect(mg,chironomid);
Matrix(1).title = 'chironomid';
Matrix(1).color = [170/255 24/255 24/255];


%% preallocate Bins
binVec=min(binEdges):binStep:max(binEdges);
binMid=binVec(1:end-1)+binStep/2;
binVec = binVec';
binMid = binMid';

plotEdges=reshape([binVec(1:end-1) binVec(2:end)]',[],1);

%% Loop though each proxy and calculate the temperature average. Include # of records contributing to average.

for i = 1:8 %step through each proxy
    
    % index of records within latitude band
    d = Matrix(i).index; %intersect(mg,banding);
    
    Matrix(i).indexRecs = d;
    
    %bin things
    
    for j=1:length(d) %step though each record in lat band
        %saving the record names in the latitude bands
        Matrix(i).names{j,1} = TS(d(j)).dataSetName;
        [bin_mean, BinTime, sem, bin_sum] =  bin_x(TS(d(j)).age,TS(d(j)).paleoData_values,binVec);
        if strcmp(TS(d(j)).interpretation1_direction,'negative') %checking the climate interpretation direction
            bin_mean = -1*bin_mean; %Flipping if interpretation direction negative
        end
        
        bin_mean=(bin_mean-nanmean(bin_mean))/nanstd(bin_mean);
        %putting the records from a band in a matrix
        Tmat(:,j) = bin_mean;
    end
    
    Tcomp = nanmean(Tmat,2);
    
    %saving all the medians for ever more
    Matrix(i).records = Tmat;
    Matrix(i).Median = Tcomp;
    Matrix(i).time = binMid;
    Matrix(i).plotEdges = plotEdges;
    
    %sample depth
    t = ~isnan(Tmat);
    t2 = sum(t,2);
    Matrix(i).sampleDepth = t2;
    
    %bootstrap sampling
    records = Matrix(i).records;
    
    p_boot = bootstrp(1000,@nanmean,records');
    
    Matrix(i).medianEnsemble = p_boot';
    
    clear Tcomp
    clear t2
    clear Tmat
 
end


%% plotting
bands = Matrix;

cc = 1.1;

style_l = {'FontName','Heveltica','FontSize',12,'FontWeight','bold'};
style_i = {'FontName','Heveltica','FontSize',10};

fig('Fig4')

for i = 1:length(bands);
    
    %get some error bars
    ci_archs(i).quant  = quantile(bands(i).medianEnsemble,[0.025 0.975],2);
    ciLo = reshape([ci_archs(i).quant(:,1) ci_archs(i).quant(:,1)]',[],1);
    ciHi = reshape([ci_archs(i).quant(:,2) ci_archs(i).quant(:,2)]',[],1);
    
    ylims       = [floor(cc*min(ciLo)), ceil(cc*max(ciHi))];
    ensembleMedian = bands(i).Median;
    
    plotScale=reshape([ensembleMedian ensembleMedian]',[],1);
    
    %start the subplotting
    subplot(round(length(bands)/2),2,i);
    
    [ax,h1,h2] = plotyy(bands(i).time,bands(i).sampleDepth,bands(i).plotEdges,plotScale,@bar,@line); hold on
    set(h1,'facecolor',rgb('Gainsboro'),'edgecolor','none')
    set(ax(1),'Ycolor',rgb('Silver'),'YAxisLocation','right','TickDir','out','YMinorTick','off')
    set(get(ax(1),'Ylabel'),'String','# records')
    
    %set sample depth plot limits and things
    set(get(ax(1),'Ylabel'),style_l{:}), set(ax(1),'yLim',[0 max(bands(i).sampleDepth)+5]), %set(ax(1),'Ytick',10:10:roundn(max(bands(i).sampleDepth),1))
    set(h2,'color',bands(i).color,'linewidth',2);
    set(ax(2),'Ycolor',bands(i).color,'YAxisLocation','left')
    set(get(ax(2),'Ylabel'),'String','z-score');
    set(get(ax(2),'Ylabel'),style_l{:}),
    set(gcf,'CurrentAxes',ax(2))
    set(ax(2),'yLim',[nanmean(plotScale)-1.5 nanmean(plotScale)+1.5])
    set(ax(2),'Ytick',round(nanmean(plotScale)-1.5):.5:round(nanmean(plotScale)+1.5))
    
    
    set(ax(1),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
    set(ax(2),'Box', 'off', 'TickDir','out','TickLength',[.02 .02], 'fontsize',11,'fontname','Heveltica')
    set(ax(1),'Xdir','reverse','Xlim',[-100 12000])
    set(ax(2),'Xdir','reverse','Xlim',[-100 12000])
    set(ax(1),'Xtick',0:2000:12000,'xminortick','on')
    
    % plot bootstrap CI
    wide = (~isnan(ciLo) & ciHi-ciLo>range(ylims)/50.0);  %
    hci = area_binFill(bands(i).plotEdges(wide)',ciLo(wide)',ciHi(wide)',bands(i).color,bands(i).color,0.1);
    
    linkaxes([ax(1),ax(2)],'x'); set(ax(1),'xLim',[-100 12000]);
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','off', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    set(ax(2),'XTick',[])
    ch = get(h1,'child'); set(ch,'EdgeAlpha',.3)
    
    %add some titles
    if i == 1
        h = title([bands(i).title ', n = ' int2str(size(bands(i).records,2)) ''],'fontname','Heveltica');
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-.02 P(3)])
        set(gcf,'color','white')
        
    else
        h = title([bands(i).title ', n = ' int2str(size(bands(i).records,2)) ' '],'fontname','heveltica');
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-0.5 P(3)])
        set(gcf,'color','white')
    end
    
    %setting some axis limits based on sample depths for looks
    if max(bands(i).sampleDepth)>=200
        set(ax(1),'Ytick',0:100:roundn(max(bands(i).sampleDepth),2))
    elseif max(bands(i).sampleDepth)>=60
        set(ax(1),'Ytick',0:20:roundn(max(bands(i).sampleDepth),1))
    else
        set(ax(1),'Ytick',0:10:roundn(max(bands(i).sampleDepth),1))
    end
    
    %laveling x-axis of the last plots in the loop
    if i == 7;
        set(get(ax(1),'Xlabel'),'String','year (BP) ')
        set(get(ax(1),'Xlabel'),style_i{:})
    end
    
    if i == 8;
        set(get(ax(1),'Xlabel'),'String','year (BP) ')
        set(get(ax(1),'Xlabel'),style_i{:})
    end
    
end

