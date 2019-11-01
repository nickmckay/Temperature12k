%% Site location map (Temperature 12k)
% Cody Routson, Oct 2019

clear, close all
load TS.mat


%% Make sure all the year columns are populated or it breaks the availability plot. Could code more elegantly. 

for i = 1:length(TS)
    TS(i).year = 1950-TS(i).age;
end

%% setting colors. 


Graph{1,1}=[8/255 46/255 114/255];  
Graph{2,1}=[221/255 84/255 46/255]; 
Graph{3,1}=[53/255 165/255 197/255];    
Graph{4,1}=[170/255 24/255 24/255]; 
Graph{5,1}= [236/255 156/255 46/255]; 
Graph{6,1}=[33/255 52/255 219/255];    
Graph{7,1}=[78/255 188/255 249/255]; 
Graph{8,1}=[128/255 54/255 168/255];         
Graph{9,1}=[50/255 127/255 81/255];  


%marker size and some text styles

style_t = {'FontName','Helvetica','FontSize',14,'FontWeight','bold'};
style_l = {'FontName','Helvetica','FontSize',12,'FontWeight','bold'};
ms = 8;


%% Index the database

ma = find(strcmp('Temp12k',{TS.paleoData_inCompilation}')); %index the compilation

mb = find(strncmp('annual',{TS.interpretation1_seasonalityGeneral}',7)); %index all annual records
mc = find(strncmp('summerOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all SummerOnly records
me = find(strncmp('winterOnly',{TS.interpretation1_seasonalityGeneral}',7)); %index all WinterOnly records

mf = [mb; mc; me]; %index of annual, summer and winter

mg = intersect(mf,ma); %intersect seasons wanted with Temp12k compilation for a USE IN index


%% Indexing proxies

aproxy = {TS.paleoData_proxyGeneral};
proxy = aproxy(mg);

%% dealing with un-assigned proxies if necessary
emptyCells = cellfun(@isempty,proxy);
proxy(emptyCells) = {'empty'};

%% put macrofossil into biophysical lump
cl = find(strcmpi('macrofossils',proxy));
proxy(cl) = {'biophysical'};

cx = find(strcmpi('microfossil',proxy));
proxy(cx) = {'other microfossil'};

%%
proxies = unique(proxy);
%% proxy and archive codes

%proxy codes
p_code = nan(size(proxy));

for i=1:length(proxies)
    id = strcmp(proxy, proxies(i));
    p_code(id) = i;
end

%archive codes. 1 through number of archives
archivee = {TS.archiveType};
archive = archivee(mg);
archives = unique(archive);
a_code = nan(size(archive));

for j=1:length(archives)
    id = strcmp(archive, archives(j));
    a_code(id) = j;
end

%% St*ff for the temporal availability plot by proxy and for archive legend
year = -10000:2000;
ny = length(year);
nvar = length(mg);
avail = NaN(ny,nvar);

for r = 1:nvar
    yearMin = round(min(TS(mg(r)).year));
    yearMax = round(max(TS(mg(r)).year));
    avail(ismember(year,[yearMin:yearMax]),r)=1;
    edgec{r} = 'k';
    
end

na = length(proxies);
nproxy = zeros(ny,na); 
pind = zeros(na,1);
for a = 1:na % loop over proxies
    nproxy(:,a) = sum(~isnan(avail(:,p_code == a)),2); %*This gets plotted in proxy temporal availability
    pind(a) = find(p_code == a,1,'first');
end

%getting archive types for legend. This can go away if we don't plot the
%archives as different shapes

nb = length(archives);
aind = zeros(nb,1);

for a = 1:nb % loop over archive types
    aind(a) = find(a_code == a,1,'first');
end

%% St*ff for temporal availability by season

% index seasons
aseason = {TS.interpretation1_seasonalityGeneral};
seasonAvail1 = aseason(ma);
emptyCell = cellfun(@isempty,seasonAvail1);
seasonAvail1(emptyCell) = {'annual'};

% Cleanup. Replace Annual with annual
axx = find(strcmpi('Annual',seasonAvail1));
seasonAvail1(axx) = {'annual'};

% replace season names with generics (winter, summer, annual)
seasonAvail = seasonAvail1;
wint1 = find(strcmpi('winterOnly',seasonAvail));
wint2 = find(strcmpi('winter+',seasonAvail));
wint3 = find(strcmpi('winter: need to combine',seasonAvail)); %these "need to combine" should be combined in newer database versions

win = [wint1, wint2, wint3];

% assigning the season lumps (renaming into lumps)
for ki = 1:length(win);
    seasonAvail{win(ki)} = 'winter';
end

sum1 = find(strcmpi('summerOnly',seasonAvail));
sum2 = find(strcmpi('summer+',seasonAvail));
sum3 = find(strcmpi('summer: need to combine',seasonAvail));

sume = [sum1, sum2, sum3];

% assigning the season lumps (renaming into lumps)
for ri = 1:length(sume);
    seasonAvail{sume(ri)} = 'summer';
end

% seasons = {'summerOnly', 'winterOnly', 'annual'};
seasons = unique(seasonAvail);

%season codes
s_code = nan(size(seasonAvail));

for i=1:length(seasons)
    id = strcmp(seasonAvail, seasons(i));
    s_code(id) = i;
    
end

nvarr = length(ma);
avails = NaN(ny,nvarr);

for r = 1:nvarr
    yearMin = round(min(TS(ma(r)).year));
    yearMax = round(max(TS(ma(r)).year));
    avails(ismember(year,[yearMin:yearMax]),r)=1;
    edgec{r} = 'k';
    
end

ns = length(seasons);
nseas = zeros(ny,ns); 
seasonal = unique(seasons);

for a = 1:ns % loop over seasons
    nseas(:,a) = sum(~isnan(avails(:,s_code == a)),2); 
    %*This gets plotted in proxy temporal availability
end


%%
%getting season types for legend
nbs = length(seasonal);
sind = zeros(nbs,1);

for ax = 1:nbs % loop over archive types
    sind(ax) = find(s_code == ax,1,'first');
end

%% lats and lons for all proxies to plot

latitude={TS.geo_latitude};
lat = (cell2mat(latitude))';

longitude={TS.geo_longitude};
lon = (cell2mat(longitude))';

p_lon = lon(mg);
p_lat = lat(mg);

%% MAP PLOT
% =============
% SET THE ST*FF
% =============

fig('Map'), clf
% orient landscape
set(gcf,'PaperPositionMode','auto')

% figure size [left bottom width height]
set(gcf, 'Position', [440   144   1000   1200])

% map pannel size
hmap=axes('Position', [.05 0.61 0.65 0.4]);

% do the map things
axesm( 'MapProjection','Robinson','MapLatLimit',[ -90 90 ],'MapLonLimit',[ -180 180 ],'MLineLocation',20,...
    'PLineLocation',20,'MeridianLabel','off','ParallelLabel','on',...
    'LabelFormat','compass','MLabelParallel','south','fontsize',12 );
tightmap;
paperscale;
axis off; framem on;

geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8])
gridm('glinestyle','.','plinelocation',20)

load coast;
linem(lat,long,'Color','k','LineWidth',2);

% =============
% PLOT LOCATIONS AND COLOR BY GROUP IN A SLOW LOOP
% =============

for r = 1:nvar
    h(r) = plotm(rand(1)*2+p_lat(r),rand(1)*2+p_lon(r),'marker','o','MarkerEdgeColor',edgec{r},'MarkerFaceColor',Graph{p_code(r),1},'linewidth',[1],'MarkerSize',[ms],'linestyle','none');

end


% =============
% PLOT LEGEND
% =============

legNames = proxies';

hl = legend(h([pind]),legNames,'location',[.79 .45 .1 .2],style_l{:});
set(hl, 'FontName', 'Palatino','box','off');


% =============
% PLOT PROXY TEMPORAL COVERAGE
% =============

%matching color map for proxies
hstack=axes('Position', [0.08 0.38 0.6 0.22]);
cmap=cell2mat(Graph(:,1));
cmap = flipud(cmap); %flipping, flopping and truncating to get the correct colors on the temporal plot to match proxy legend
cmap = cmap(end-size(nproxy,2)+1:end,:);

colormap(gca,cmap);

nproxy = fliplr(nproxy);
area(1950-year,nproxy,'EdgeColor','w'), set(gca,'YAxisLocation','Right');
set(gca, 'XDir', 'reverse')
xlim([-60 12000])
fancyplot_deco('','year (BP)','# records');
title('Temporal Availability',style_t{:})

set(gca,'xtick',0:500:12000)

set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','off','YMinorTick','on', 'YGrid','on')


% =============
%  PLOT SEASON TEMPORAL COVERAGE
% =============

%[left bottom width height]
hstack=axes('Position', [0.08 0.07 0.6 0.21]);

%season color map
col(1,1:3) = rgb('black');
col(2,1:3) = rgb('firebrick')';
col(3,1:3) = rgb('steelblue')';

colormap(gca, col)

HD = area(1950-year,nseas,'EdgeColor','w');
set(gca,'YAxisLocation','Right');
set(gca, 'XDir', 'reverse')
xlim([-60 12000])
fancyplot_deco('','year (BP)','# records');

set(gca,'xtick',0:500:12000)


set(gca,'Box', 'off', 'TickDir','out','TickLength',[.02 .02],'XMinorTick','off','YMinorTick','on', 'YGrid','on')
title('Seasonal Availability',style_t{:})


hx = legend(HD,seasonal','location',[.79 .1 .1 .2],style_l{:});
set(hx, 'FontName', 'Palatino','box','off');

