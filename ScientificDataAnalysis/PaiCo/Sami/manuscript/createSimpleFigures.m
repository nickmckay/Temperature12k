function createSimpleFigures

mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'results' filesep 'simpleexperiment_new.mat']);
mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'manuscript_colors.mat']);

names = fieldnames(results{1}{1});
times = 1:500;
calibLength = 100;

figure(1); clf; hold all;

means = results{1}{1};
jobResult = cell2mat(results{1});
for j = 1:numel(names)
    means.(names{j}) = mean(cat(1,jobResult.(names{j})),1);
end
methods = {'target' 'PaiCo' 'MOM' 'OLS' 'PCReg' 'RegEM' 'LNA' 'BARCAST'};
methodColors = [[1 0 0 ];darkblue;green;cyan;brown;magenta;lightgray;lightblue];

for i = 1:numel(names)
    plot(times, means.(names{i}),'color',methodColors(strcmpi(methods,names{i}),:));
end
legend(names,'location','southwest');
axis tight;
axs = axis;
line(times(end-calibLength+1)*[1 1], axs(3:4),'color','black');

squarepage([8,8]);
print('-dpdf',[mpath filesep 'figures' filesep 'simple_result_base.pdf']);

comparisonMask = times <= 100;



for i = 2:5
    names = fieldnames(results{i}{1});
    jobResult = cell2mat(results{i});
    errors = nan(size(jobResult,1),size(jobResult,2),numel(names));
    for k = 1:size(jobResult,1)
        for j = 1:numel(names)
            singleResult = jobResult(k,:);
            res = cat(1,singleResult.(names{j}));
            errors(k,:,j) = mean(sqrt(mean((res(:,comparisonMask)-jobs(i).startHeight).^2))/jobs(i).startHeight);
        end
    end
    
    err = squeeze(mean(errors,2));
    
    x = jobs(i).param;
    if (iscell(x))
        x = x{:};
    end
    if (i == 4)
        figure((i-1)*2); clf; hold all;
        for j = 1:size(jobs(i).noisevar,1)
            plot(jobs(i).noisevar(j,:));
        end        
        axis tight;
        squarepage([8,8]);        
        print('-dpdf',[mpath filesep 'figures' filesep 'simple_job_' num2str(i) '.pdf']);                
    elseif (i == 5)
        figure((i-1)*2); clf; hold all;        
        sr = jobResult(1,1);
        y = linspace(min(sr.target), max(sr.target),1e2);
        mask = (y >= -1) & (y <= 1);
        for j = 1:size(jobs(i).transfer,1)
            t = jobs(i).transfer{j}(y);
            plot(y,t/(max(t(mask))-min(t(mask))));
        end        
        axis tight;
        squarepage([8,8]);        
        print('-dpdf',[mpath filesep 'figures' filesep 'simple_job_' num2str(i) '.pdf']);        
    end
    
    figure((i-1)*2+1); clf; hold all;
    for j = 1:numel(names)
        plot(x,err(:,j), 'color', methodColors(strcmpi(names{j}, methods),:));
    end
    axis tight;
    axs = axis;
    
    if (i == 2)
        axs(4) = 0.45;
    elseif (i == 3)
        axs(4) = .28;
    end
    axis(axs);
    
    squarepage([8,8]);
    print('-dpdf',[mpath filesep 'figures' filesep 'simple_result_' num2str(i) '.pdf']);
end


