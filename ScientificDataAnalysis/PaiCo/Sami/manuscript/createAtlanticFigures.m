function createAtlanticFigures

mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'manuscript_colors.mat']);
load([mpath filesep 'results' filesep 'arcticatlantic.mat']);
mpath = fileparts(mfilename('fullpath'));

if (~exist('m_proj'))
    error('Requires the MMap package for Matlab. http://www2.ocgy.ubc.ca/~rich/map.html ... Edit CD: Link outdated. New link: https://www.eoas.ubc.ca/~rich/map.html');
end

coord = zeros(numel(records),2);
typeMap = java.util.HashMap;
typeMap.put('Ice core',0);
typeMap.put('Marine sediment',1);
typeMap.put('Historic',2);
typeMap.put('Speleothem',3);
typeMap.put('Tree ring',4);
typeMap.put('Lake sediment',5);
types = nan(numel(records),1);
for i = 1:numel(records)
    coord(i,1) = records{i}.latitude;
    coord(i,2) = records{i}.longitude;  
    t = typeMap.get(records{i}.proxytype);
    if (isempty(t))
        disp(['new type: ' records{i}.proxytype ' ' num2str(typeMap.size()+1)]);         
        typeMap.put(records{i}.proxytype,typeMap.size());
        t = typeMap.get(records{i}.proxytype);
    end
    types(i) = t;
end


%% Arctic map of record locations
figure(1); clf; hold all;
m_proj('lambert','long',[-60 40],'lat',[55 90]);
m_coast('patch',[1 .85 .7]);

% Sort records
uq = unique(types);
newind = zeros(numel(types),1);
pos = 0;
for i = 1:numel(uq)
    tind = find(types == uq(i));
    [x,y] = m_ll2xy(coord(tind,2), coord(tind,1));
    [~,ind] = sort(y,'ascend');
    newind((1:numel(ind))+pos) = tind(ind);
    pos = pos + numel(ind);
end
ind = newind;
coord = coord(ind,:);
data.proxy = data.proxy(ind);
records = records(ind);
types = types(ind);

m_grid('xtick',5,'tickdir','out','ytick',[60 70 80],'linest','-');
set(gca,'FontSize',8);
set(findobj(gca,'type','text'),'fontsize',10);
m_plot([linspace(-50,30,100) 30 -50], [60*ones(1,100) 90 60],'k-', 'linewidth',1);
typeMarkers = '>+sop^';

for i = 1:numel(uq)    
    mask = types == uq(i);
    m_plot(coord(mask,2),coord(mask,1),['k' typeMarkers(i)],'markersize',4,'markerfacecolor',[0,0,0]);
    
end
labels = cellfun(@(x)(num2str(x)), num2cell(1:numel(records)),'UniformOutput',false);
m_text(coord(:,2), coord(:,1), labels);
squarepage([10 10]);
print('-dpdf',[mpath filesep 'figures' filesep 'arctic_proxies.pdf']);

%% Temporal extent of proxies
years = 0:2000;
colors = rand(numel(uq),3);

figure(2); clf; hold all;
margin = 0.0;
handles = nan(1,numel(uq));
for i = 1:numel(records)
    record = records{i};
    colorind = types(i)+1;
    if (size(record.content,2)==4)
        Y = [0 0 1 1]*(1-margin)+margin/2-i;
        for j = 1:size(record.content,1)
            cover = years(years >= record.content(j,3) & years < record.content(j,4));
            if (isempty(cover))
                continue;
            end
            xmin = min(cover);
            xmax = max(cover)+1;
            h = patch([xmin xmax xmax xmin], Y, colors(colorind,:));
        end
        if (isnan(handles(colorind)))
            handles(colorind) = h;
        end
    else
        Y = [0 0 1 1]*(1-margin)+margin/2-i;
        cover = record.content(:,1)';
        mask  = [false ismember(cover,years) false];
        indinc = find(diff(mask) == 1);
        inddec = find(diff(mask) == -1)-1;
        for j = 1:numel(indinc)
            xmin = min(cover(indinc(j)));
            xmax = max(cover(inddec(j)))+1;
            h = patch([xmin xmax xmax xmin], Y, colors(colorind,:));
        end
        if (isnan(handles(colorind)))
            handles(colorind) = h;
        end
        
    end
end
set(gca,'ytick',[]);
xlabel('Year AD');
axis tight;
print('-dpdf',[mpath filesep 'figures' filesep 'temporalextent.pdf']);

%% Plot proxy information in tabular for to be inserted in manuscript
for i = 1:numel(records)
    fprintf('%d\t%s\t%s\t%s\t%5.2f\t%5.2f\t%s\t%d--%d\t%s\t%5.2f\t%s\n',i, records{i}.region, records{i}.area, records{i}.site, coord(i,1), coord(i,2), records{i}.proxytype, min(records{i}.content(:,1)), max(records{i}.content(:,1)), records{i}.measurement, records{i}.resolution, records{i}.id);
end

%% Plot original results with focus to 1800-2000
figure(3); clf;
hold all;
msk = result.times >= 1800;
errorhandle = area(result.times(msk), [result.signal(msk)'-result.noisestd 2*ones(nnz(msk),1)*result.noisestd],'EdgeColor',errorColor);
set(errorhandle(1),'FaceColor',[1 1 1]);
set(errorhandle(2),'FaceColor',errorFaceColor);
set(get(get(errorhandle(1),'Annotation'),'LegendInformation'), 'IconDisplayStyle','off');

if (isfield(result,'resample'))
    resample = sort(result.resample.signals,1);
    fpmin = resample(round(size(resample,1)*0.05),:);
    fpmax = resample(round(size(resample,1)*0.95),:);    
    fpmean = resample(round(size(resample,1)/2),:);
    handles(4) = plot(result.times, fpmean, '--', 'color', magenta,'linewidth',1);
    handles(5) = plot(result.times, fpmin, ':', 'color', magenta,'linewidth',1);
    plot(result.times, fpmax, ':', 'color', magenta,'linewidth',1);
end

instru = nanmean(data.instrumental.data,1);

paicohandle = plot(result.times(msk), result.signal(msk),'-','color',darkblue);
msk = data.instrumental.times >= 1800;
instruhandle = plot(data.instrumental.times(msk), instru(msk),'-','color',red);
legend([instruhandle paicohandle errorhandle(2)], {'Instrumental' 'PaiCo' 'Estimated error'},'location','northwest');
axis tight;
axs = axis;
axis([1800, axs(2:end)]);
    
xlabel('Years (AD)');
ylabel('Anomaly (Celsius)');

squarepage([20 5]);
print('-dpdf',[mpath filesep 'figures' filesep 'paico_focus.pdf']);

%% Plot smoothed results with different error estimates
figure(4); clf;
span = 10;

linewidth = 0.5;
hold all;
handles = [];

[times, signal] = aggregate(result.times, result.signal, span);

handles(3) = patch([times(:);flipud(times(:))], [signal(:)-result.noisestd/span;flipud(signal(:))+result.noisestd/span],1);
set(handles(3),'EdgeColor',errorColor,'FaceColor',errorFaceColor);

if (isfield(result,'resample'))
    resample = sort(result.resample.signals,1);
    [times, resample] = aggregate(result.times, resample, span);
    fpmin = resample(round(size(resample,1)*0.05),:);
    fpmax = resample(round(size(resample,1)*0.95),:);    
    fpmean = resample(round(size(resample,1)/2),:);
    handles(4) = plot(times, fpmean, '--', 'color', magenta,'linewidth',linewidth);
    handles(5) = plot(times, fpmin, ':', 'color', magenta,'linewidth',linewidth);
    plot(times, fpmax, ':', 'color', magenta,'linewidth',linewidth);
end

instru = nanmean(data.instrumental.data,1);
[times, signal] = aggregate(data.instrumental.times, instru, span);
handles(1) = plot(times, signal,'-','color',red,'linewidth',linewidth);
[times, signal] = aggregate(result.times, result.signal, span);
handles(2) = plot(times, signal,'-','color', darkblue,'linewidth',linewidth);
lh = legend(handles, {'Instrumental' 'PaiCo' 'Estimated error' 'Bootstrap median' 'Bootstrap 90%'},'location','northoutside','orientation','horizontal');
plot(result.times, result.signal);
axis tight;

xlabel('Years (AD)');
ylabel('Anomaly (Celsius)');

squarepage([20 5]);

print('-dpdf',[mpath filesep 'figures' filesep 'paico_error_comparison.pdf']);


%% Plot aggregated decadal coparison to previous studies
figure(5); clf;
hold all;


handles = [];
pdata = xlsread('kaufman.xlsx');
ptimes = pdata(:,1);
pdata = pdata(:,2:end)';
% [times, pdata] = aggregate(result.times, pdata, span);
colors = [cyan; green; lightgray; brown;];
for i = 1:size(pdata,1)
    handles(i+2) = plot(ptimes, pdata(i,:), '-', 'color', colors(i,:));
end
    
linewidth = 1;
instru = nanmean(data.instrumental.data,1);
[times, signal] = aggregate(data.instrumental.times, instru, span);
handles(1) = plot(times, signal,'-','color',red,'linewidth',linewidth);
[times, signal] = aggregate(result.times, result.signal, span);
handles(2) = plot(times, signal,'-','color', darkblue,'linewidth',linewidth);
lh = legend(handles, {'Instrumental' 'PaiCo' 'Mann et al. 2008 (CPS)' 'Mann et al. 2008 (EIV)' 'Moberg et al. 2005' 'Kaufmann et al. 2009'},'location','north','orientation','horizontal','box','off');

xlabel('Years (AD)');
ylabel('Anomaly (Celsius)');
axis tight;


squarepage([20 5]);
print('-dpdf',[mpath filesep 'figures' filesep 'paico_reconst_comparison.pdf']);


% Calculate the most pronounced hot and cool episodes

wind = 10:10:400;
d = zeros(numel(wind),numel(result.signal));
hotcool = zeros(1,numel(result.signal));
for i = 1:numel(wind)
    window = 1:wind(i);
    for j = 1:(numel(result.signal)-wind(i)+1);
        X = [ones(wind(i),1) result.times(j-1+window)'];
        alpha = X\result.signal(j-1+window)';
        d(i,j) = alpha(2);
    end    
    
    [~,ind] = sort(abs(d(i,:)),'descend');
    
    for j = 1:3
        X = [ones(numel(window),1) result.times(ind(j)-1+window)'];
        alpha = X\result.signal(ind(j)-1+window)';
        hotcool(ind(j)-1+window) = hotcool(ind(j)-1+window) + sign(alpha(2));
    end   
end
figure(6); clf; hold all;
threshold = 1;
mask = abs(hotcool) > threshold;
hotcool(mask) = threshold*sign(hotcool(mask));
minmax = [min(result.signal) max(result.signal)];
i = 1;
hotcools = zeros(numel(hotcool),2);
pos = 0;
while (i < numel(hotcool))
    j = find(hotcool((i+1):end) ~= hotcool(i),1,'first');
    if (isempty(j))
        if (all(hotcool((i+1):end) == hotcool(i)))
            j = numel(hotcool)-i;
        else
            j = 1;
        end
    end
    a = patch([result.times(i) result.times(i+j) result.times(i+j) result.times(i)],...
        [minmax(1) minmax(1) minmax(2) minmax(2)], hotcool(i));
    
    if (hotcool(i) ~= 0)
        pos = pos + 1;
        hotcools(pos,:) = [i (i+j)];
        alpha = [ones(j+1,1) (0:j)']\result.signal(i:(i+j))';
        plot([result.times(i) result.times(i+j)], alpha(1)+[0 j*alpha(2)],'g-','linewidth',1);
    end
    set(a,'EdgeColor','none');
    
    i = i+j;
end

% Save positions of coolings and warmings for uncertainty analysis
hotcools = hotcools(1:pos,:);
save([mpath filesep 'results' filesep 'hotcools.mat'],'hotcools');

handles(1) = plot(result.times, result.signal,'k');
[stime, ssig] = aggregate(result.times, result.signal, 10);
handles(2) = plot(stime, ssig, 'm','linewidth',2);
colormap([linspace(0,1,100) ones(1,100); linspace(0,1,100) linspace(1,0,100); ones(1,100) linspace(1,0,100)]');


xlabel('Years (AD)');
ylabel('Anomaly (Celsius)');
legend(handles,{'PaiCo','PaiCo, 10-yr bins'});
axis tight;


squarepage([20 5]);
print('-dpdf',[mpath filesep 'figures' filesep 'paico_hotcool.pdf']);

%% Statistical analysis of signals

signalR = result.signal(ismember(result.times, data.instrumental.times));
signalI = data.instrumental.data(ismember(data.instrumental.times, result.times));
signalR = signalR(:);
signalI = signalI(:);
times = result.times(ismember(result.times, data.instrumental.times));

[~, signalRd] = aggregate(times, signalR, 10);
[~, signalId] = aggregate(times, signalI, 10);
signalRd = signalRd(:);
signalId = signalId(:);

p = polydata(ar(signalR, 1));

scorr = zeros(1e4, 2);
for i = 1:size(scorr,1)    
    t = randn(size(signalR));
    t(2:end) = t(2:end) - p(2)*t(1:(end-1));
    
    scorr(i,1) = corr(signalI, t);
    [~, t] = aggregate(times, t, 10);
    t = t(:);
    scorr(i,2) = corr(signalId, t); 
end
c = corr(signalR, signalI);
disp(['Annual correlation with instrumental = ' num2str(c) ' p=' num2str((sum(scorr(:,1) >= c)+1)/size(scorr,1))]);
cd = corr(signalRd, signalId);
disp(['Decadal correlation with instrumental = ' num2str(cd) ' p=' num2str((sum(scorr(:,2) >= cd)+1)/size(scorr,1))]);

X = [ones(numel(result.times),1) result.times(:)];
alpha = X\result.signal(:);
disp(['Coefficient of linear fit to PaiCo = ' num2str(1e3*alpha(2)) ' C/1000yrs']);
