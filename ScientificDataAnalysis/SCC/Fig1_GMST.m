clear
close all
%%
load sampleDepth.mat

%%
cd ./final_data

files = dir('*.csv');

latName{1} = '60°N to 90°N';
latName{2} = '30°N to 60°N';
latName{3} = '0°N to 30°N';
latName{4} = '30°S to 0°S';
latName{5} = '60°S to 30°S';
latName{6} = '90°S to 60°S';

%% Load the csv files for the latitudinal recontructions
names = {files(:).name};
names = names';

IndexC = strfind(names,'CPS_latband');
ind = find(not(cellfun('isempty',IndexC)));

delimiterIn = ',';
headerlinesIn = 0;

for i = 1:6
    CPS(i).title = files(ind(i)).name;
    B = importdata(files(ind(i)).name,delimiterIn,headerlinesIn);
    CPS(i).ensemble =  B(:,2:end);
    CPS(i).time = B(:,1);
end


% load DCC
clear IndexC
clear ind
clear B
clear i

IndexC = strfind(names,'DCC_latband');
ind = find(not(cellfun('isempty',IndexC)));

delimiterIn = ',';
headerlinesIn = 0;

for i = 1:6
    DCC(i).title = files(ind(i)).name;
    B = importdata(files(ind(i)).name,delimiterIn,headerlinesIn);
    DCC(i).ensemble =  B(:,2:end);
    DCC(i).time = B(:,1);
end


% load GAM
clear IndexC
clear ind
clear B
clear i

IndexC = strfind(names,'GAM_latband');
ind = find(not(cellfun('isempty',IndexC)));

delimiterIn = ',';
headerlinesIn = 0;

for i = 1:6
    GAM(i).title = files(ind(i)).name;
    B = importdata(files(ind(i)).name,delimiterIn,headerlinesIn);
    GAM(i).ensemble =  B(:,2:end);
    GAM(i).time = B(:,1);
end


% load PAI
clear IndexC
clear ind
clear B
clear i

IndexC = strfind(names,'PAI_latband');
ind = find(not(cellfun('isempty',IndexC)));

delimiterIn = ',';
headerlinesIn = 0;

for i = 1:6
    PAI(i).title = files(ind(i)).name;
    B = importdata(files(ind(i)).name,delimiterIn,headerlinesIn);
    PAI(i).ensemble =  B(:,2:end);
    PAI(i).time = B(:,1);
end

% load SCC
clear IndexC
clear ind
clear B
clear i

IndexC = strfind(names,'SCC_latband');
ind = find(not(cellfun('isempty',IndexC)));

delimiterIn = ',';
headerlinesIn = 0;
%%
for i = 1:6
    SCC(i).title = files(ind(i)).name;
    B = importdata(files(ind(i)).name,delimiterIn,headerlinesIn);
    SCC(i).ensemble =  B(:,2:end);
    SCC(i).time = B(:,1);
    
    clear B
end
%%
cd ..



%%
plotTime = reshape([SCC(1).time SCC(1).time]',[],1);

figure(1)
clf

a1 = 1:6:37;

style_l = {'FontName','Helvetica','FontSize',12};
cols = brewermap(6,'dark2');

set(gcf,'color','white')

%============== Composite No Scale (SCC)

for i = 1:6;
    
    ensembleMedian = nanmedian(SCC(i).ensemble,2);
    
    plotScale=reshape([ensembleMedian ensembleMedian]',[],1);
    
    ci_archs.quant  = quantile(SCC(i).ensemble,[0.025 0.975],2);
    ciLo = ci_archs.quant(:,1);
    ciHi = ci_archs.quant(:,2);
    
    ciHi = reshape([ciHi ciHi]',[],1);
    ciLo = reshape([ciLo ciLo]',[],1);
    
    subplot(6,6,a1(i));
    
    h1 = plot(plotTime,plotScale); hold on
    set(h1,'color',cols(i,:),'linewidth',2);
    set(gca,'Ycolor',cols(i,:),'YAxisLocation','left')
    set(gcf,'CurrentAxes',gca)
    set(gca,'yLim',[-4 4])
    set(gca,'Ytick',-4:1:4)
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    set(gca,'Xdir','reverse','Xlim',[-100 12000])
    set(gca,'Xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000])
    
    hci = area_binFill(plotTime',ciLo',ciHi',cols(i,:),cols(i,:));
    
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','off','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    
    h = text(10000,-2,[latName(i)],'fontweight','normal');
    
    if i<6
        set(gca,'xticklabel',[])
    end
    
    if i == 1
        h = title('SCC','fontweight','normal');
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-.5 P(3)])
        
    else if i == 3
            set(get(gca,'Ylabel'),'String','temperature (°C)')
            set(get(gca,'Ylabel'),style_l{:})
            
        end
    end
    
end


a2 = 2:6:38;

%============== DCC

for i = 1:6;
    
    ensembleMedian = nanmedian(DCC(i).ensemble,2);
    
    plotScale=reshape([ensembleMedian ensembleMedian]',[],1);
    
    ci_archs.quant  = quantile(DCC(i).ensemble,[0.025 0.975],2);
    ciLo = ci_archs.quant(:,1);
    ciHi = ci_archs.quant(:,2);
    
    ciHi = reshape([ciHi ciHi]',[],1);
    ciLo = reshape([ciLo ciLo]',[],1);
    
    subplot(6,6,a2(i));
    
    h1 = plot(plotTime,plotScale); hold on
    set(h1,'color',cols(i,:),'linewidth',2);
    set(gca,'Ycolor',cols(i,:),'YAxisLocation','left')
    set(gcf,'CurrentAxes',gca)
    set(gca,'yLim',[-4 4])
    set(gca,'Ytick',-4:1:4)
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    set(gca,'Xdir','reverse','Xlim',[-100 12000])
    set(gca,'Xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000])
    
    hci = area_binFill(plotTime',ciLo',ciHi',cols(i,:),cols(i,:));
    
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','off','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    
    set(gca,'Yticklabel',[])
    if i<6
        set(gca,'xticklabel',[])
    end
    
    if i == 1
        h = title('DCC','fontweight','normal');
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-.5 P(3)])
    end
    
end



a4 = 3:6:39;


%============== GAM

for i = 1:6;
    
    ensembleMedian = nanmedian(GAM(i).ensemble,2);
    
    plotScale=reshape([ensembleMedian ensembleMedian]',[],1);
    
    ci_archs.quant  = quantile(GAM(i).ensemble,[0.025 0.975],2);
    ciLo = ci_archs.quant(:,1);
    ciHi = ci_archs.quant(:,2);
    
    ciHi = reshape([ciHi ciHi]',[],1);
    ciLo = reshape([ciLo ciLo]',[],1);
    
    subplot(6,6,a4(i));
    
    h1 = plot(plotTime,plotScale); hold on
    set(h1,'color',cols(i,:),'linewidth',2);
    set(gca,'Ycolor',cols(i,:),'YAxisLocation','left')
    set(gcf,'CurrentAxes',gca)
    set(gca,'yLim',[-4 4])
    set(gca,'Ytick',-4:1:4)
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    set(gca,'Xdir','reverse','Xlim',[-100 12000])
    set(gca,'Xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000])
    
    hci = area_binFill(plotTime',ciLo',ciHi',cols(i,:),cols(i,:));
    
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','off','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    set(gca,'Yticklabel',[])
    if i<6
        set(gca,'xticklabel',[])
    end
    
    if i == 1
        h = title('GAM','fontweight','normal');
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-.5 P(3)])
    end
    
end


%============== paiCo
a5 = 4:6:40;
for i = 1:6;
    
    ensembleMedian = nanmedian(PAI(i).ensemble,2);
    
    plotScale=reshape([ensembleMedian ensembleMedian]',[],1);
    
    ci_archs.quant  = quantile(PAI(i).ensemble,[0.025 0.975],2);
    ciLo = ci_archs.quant(:,1);
    ciHi = ci_archs.quant(:,2);
    
    ciHi = reshape([ciHi ciHi]',[],1);
    ciLo = reshape([ciLo ciLo]',[],1);
    
    subplot(6,6,a5(i));
    
    h1 = plot(plotTime,plotScale); hold on
    set(h1,'color',cols(i,:),'linewidth',2);
    set(gca,'Ycolor',cols(i,:),'YAxisLocation','left')
    set(gcf,'CurrentAxes',gca)
    set(gca,'yLim',[-4 4])
    set(gca,'Ytick',-4:1:4)
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    set(gca,'Xdir','reverse','Xlim',[-100 12000])
    set(gca,'Xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000])
    
    hci = area_binFill(plotTime',ciLo',ciHi',cols(i,:),cols(i,:));
    
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','off','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    
    set(gca,'Yticklabel',[])
    
    if i<6
        set(gca,'xticklabel',[])
    end
    
    if i == 1
        h = title('PAI','fontweight','normal');
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-.5 P(3)])
    end
end



%===================CPS
a6 = 5:6:41;

for i = 1:6;
    
    ensembleMedian = nanmedian(CPS(i).ensemble,2);
    
    plotScale=reshape([ensembleMedian ensembleMedian]',[],1);
    
    ci_archs.quant  = quantile(CPS(i).ensemble,[0.025 0.975],2);
    ciLo = ci_archs.quant(:,1);
    ciHi = ci_archs.quant(:,2);
    
    ciHi = reshape([ciHi ciHi]',[],1);
    ciLo = reshape([ciLo ciLo]',[],1);
    
    subplot(6,6,a6(i));
    
    h1 = plot(plotTime,plotScale); hold on
    set(h1,'color',cols(i,:),'linewidth',2);
    set(gca,'Ycolor',cols(i,:),'YAxisLocation','left')
    set(gcf,'CurrentAxes',gca)
    set(gca,'yLim',[-4 4])
    set(gca,'Ytick',-4:1:4)
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    set(gca,'Xdir','reverse','Xlim',[-100 12000])
    set(gca,'Xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000])
    
    hci = area_binFill(plotTime',ciLo',ciHi',cols(i,:),cols(i,:));
    
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','off','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    
    set(gca,'Yticklabel',[])
    
    if i<6
        set(gca,'xticklabel',[])
    end
    
    if i == 1
        h = title('CPS','fontweight','normal');
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-.5 P(3)])
    end
end


%=================== Sample Depth
a7 = 6:6:42;


for i = 1:6;
    
    subplot(6,6,a7(i));
    
    h1 = bar(sampleDepth(i).time,sampleDepth(i).allSampleDepth); hold on
    h2 = bar(sampleDepth(i).time,sampleDepth(i).calSampleDepth); hold on
    
    set(h1,'facecolor',[.5 .5 .5],'edgecolor',[.5 .5 .5])
    set(h2,'facecolor',[.8 .8 .8],'edgecolor','none')
    set(gca,'Ycolor',rgb('Silver'),'YAxisLocation','right','TickDir','out','YMinorTick','on')
    
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02])
    set(gca,'Xdir','reverse','Xlim',[-100 12000])
    set(gca,'Xtick',[0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000])
    
    set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','on','YMinorTick','on', 'YGrid','on')
    set(gca,'XColor' , [.3 .3 .3], 'LineWidth', 1);
    ch = get(h1,'child'); set(ch,'EdgeAlpha',.3)
    cp = get(h2,'child'); set(cp,'EdgeAlpha',.3)
    
    if i == 1
        set(gca,'yLim',[0 150])
        set(gca,'Ytick',[50 100 150])
        set(gca,'xticklabel',[])
        h = title('Sample Depth','fontweight','normal');
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-.5 P(3)])
        
    else if i == 2
            set(gca,'yLim',[0 250])
            set(gca,'Ytick',[50 100 150 200 250])
            set(gca,'xticklabel',[])
        else if i == 3
                set(gca,'yLim',[0 70])
                set(gca,'Ytick',[20 40 60])
                set(gca,'xticklabel',[])
                set(get(gca,'Ylabel'),'String','# records')
                set(get(gca,'Ylabel'),style_l{:})
                
            else if i == 4
                    set(gca,'yLim',[0 70])
                    set(gca,'Ytick',[20 40 60])
                    set(gca,'xticklabel',[])
                else if i == 5
                        set(gca,'yLim',[0 35])
                        set(gca,'Ytick',[10 20 30])
                        set(gca,'xticklabel',[])
                    else if i == 6
                            set(gca,'yLim',[0 20])
                            set(gca,'Ytick',[5 10 15 20])
                            
                        end
                    end
                end
            end
        end
    end
end






