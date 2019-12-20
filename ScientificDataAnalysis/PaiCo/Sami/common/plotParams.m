function plotParams(params)

clf;
fields = fieldnames(params);
num = 0;
for i = 1:numel(fields)
    if (iscell(params(1).(fields{i})))
        continue;
    end    
    num = num + numel(params(1).(fields{i}));
end

n = floor(sqrt(num));
m = ceil(num/n);
index = 0;
lgnd = cell(num,1);
hold all;
for i = 1:numel(fields)
    if (iscell(params(1).(fields{i})))% || strcmp(fields{i},'sigma2_Z')
        continue;
    end
    data = cell2mat({params(:).(fields{i})});    
    subplot(numel(fields),1,i);
    hold all;
    for j = 1:size(data,1)
%         subplot(n,m,index);       
%         hold all;        
%         plot(data(j,:));

        plot(data(j,:));
        index = index + 1;
        lgnd{index} = [fields{i} ' ' num2str(j)];

%         
%         title([fields{i} ' ' num2str(j)]);
    end
    text(1,max(max(data)),fields{i});
    axis tight;
end

% legend(lgnd(1:index),'location','EastOutside');