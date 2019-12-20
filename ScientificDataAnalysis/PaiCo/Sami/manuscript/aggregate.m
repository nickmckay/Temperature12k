function [times, data] = aggregate(times, data, span)

if (size(data,1) == numel(data))
    data = data';
end

% Binned average
[uq,~,loc] = unique(floor(times/span));
T = zeros(size(data,2),numel(uq));
T(sub2ind(size(T),1:numel(loc), loc')) = 1;
T = T./repmat(sum(T,1),size(T,1),1);
data = data*T;
times = times*T;

data = data(:,1:numel(uq));
times = times(1:numel(uq));