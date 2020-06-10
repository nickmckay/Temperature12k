function [bin_mean, BinTime, sem, bin_sum] =  bin_x(Xtime,X,bins)
%[bin_mean, BinTime, standard error of the mean] =  bin_x(Xtime,X,bins);


binstep = bins(2)-bins(1);

for i = 1:length(bins)-1
    q = find(Xtime >= bins(i) & Xtime<bins(i+1));
    
    bin_mean(i,:) = [nanmean(X(q,:),1)];
    bin_sum(i,:)  = [nansum(X(q,:),1)];
    BinTime(i,1)= nanmean(Xtime(q));
    if length(find(~isnan(X(q,:)))) == 1
        sem(i,:)=NaN;
    else
    sem(i,:) = std(X(q,:))./sqrt(length(find(~isnan(X(q,:)))));
    end
end
end
