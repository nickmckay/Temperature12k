load('C:\Users\daetwyler\Desktop\Climate12k\Matlab_related\PaiCo_12k\Results\Annual_Bin200_TargetCPS_Ref1911to1995.mat')
load('C:\Users\daetwyler\Desktop\Climate12k\Matlab_related\PaiCo_12k\Results\Summer_Bin200_TargetCPS_Ref1911to1995.mat')
load('C:\Users\daetwyler\Desktop\Climate12k\Matlab_related\PaiCo_12k\Results\Annual_Bin150_TargetCPS_Ref1911to1995.mat')
load('C:\Users\daetwyler\Desktop\Climate12k\Matlab_related\PaiCo_12k\Results\Summer_Bin150_TargetCPS_Ref1911to1995.mat')
load('C:\Users\daetwyler\Desktop\Climate12k\Matlab_related\PaiCo_12k\Results\Annual_Bin100_TargetCPS_Ref1911to1995.mat')
load('C:\Users\daetwyler\Desktop\Climate12k\Matlab_related\PaiCo_12k\Results\Summer_Bin100_TargetCPS_Ref1911to1995.mat')

figure(1)
plot(PaiCoAnnual(1).times,PaiCoAnnual(1).signal)
hold on 
for i = 2:numel(PaiCoAnnual)
    plot(PaiCoAnnual(i).times,PaiCoAnnual(i).signal)
end
legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S'}, 'Location', 'best')
hold off

figure(2)
plot(PaiCoSummer(1).times,PaiCoSummer(1).signal)
hold on 
for i = 2:numel(PaiCoSummer)
    plot(PaiCoSummer(i).times,PaiCoSummer(i).signal)
end
legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S'}, 'Location', 'best')
hold off

%% with smoothing
figure(5)
plot(PaiCoSummer(1).times,smooth(PaiCoSummer(1).signal,100,'loess'))
hold on 
for i = 2:numel(PaiCoSummer)
    plot(PaiCoSummer(i).times,smooth(PaiCoSummer(i).signal,100,'loess'))
end
legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S'}, 'Location', 'best')
hold off

figure(6)
plot(PaiCoAnnual(1).times,smooth(PaiCoAnnual(1).signal,100,'loess'))
hold on 
for i = 2:numel(PaiCoAnnual)
    plot(PaiCoAnnual(i).times,smooth(PaiCoAnnual(i).signal,100,'loess'))
end
legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S'}, 'Location', 'best')
hold off


%% Some testplotting

%% Just get data...
[test1, test2, test3] = PaiCoDataByLat(TSG_annual, nLatBands, doBinning, binWidth, instrTarget);
[test1, test2, test3] = PaiCoDataByLat(TSG_summer, nLatBands, doBinning, binWidth, instrTarget);

figure(3)
plot(test2(1).age*(-1)+1950,test2(1).values)
hold on 
for i = 2:numel(test2)
    plot(test2(i).age*(-1)+1950,test2(i).values)
end
hold off

%test3(1) corresponds to the northern most lat. band
figure(4)
plot(test3(1).proxy{1,1}.times,test3(1).proxy{1,1}.data)
hold on 
for i = 2:numel(test3(1).proxy)
    plot(test3(1).proxy{1,i}.times,test3(1).proxy{1,i}.data)
end
hold off

TSG_annualtest = TSG_annual(find(strncmpi('degC',{TSG_annual.paleoData_units}',1))); %index all calibrated
[test1, test2, test3] = PaiCoDataByLat(TSG_annualtest, nLatBands, doBinning, binWidth, instrTarget);
test4 = test3;
testtime = -9900:binWidth:1900;
testmedian = zeros(6,numel(testtime));
for i = 1:6
    for j = 1:numel(test4(i).proxy)
        %subtract mean from each proxy records
        test4(i).proxy{1,j}.data = test3(i).proxy{1,j}.data - mean(test3(i).proxy{1,j}.data);
    end
    testmat = zeros(numel(test4(i).proxy),numel(testtime));
    for k = 1:numel(testtime)
        for l = 1:numel(test4(i).proxy)
            if isempty(find(test4(i).proxy{1,l}.times == testtime(k)))
                testmat(l,k) = NaN;
            else
                testmat(l,k) = test4(i).proxy{1,l}.data(find(test4(i).proxy{1,l}.times == testtime(k)));
            end
        end
    end
    testmedian(i,:) = median(testmat,1,'omitnan');
end
figure(5)
plot(testtime,testmedian(1,:))
hold on 
for i = 2:6
    plot(testtime,testmedian(i,:))
end
legend({'60-90°N', '30-60°N', '0-30°N', '0-30°S', '30-60°S', '60-90°S'}, 'Location', 'best')
hold off




figure(5)
plot(test2(1).age*(-1)+1950,PERCENTRANK(test2(1).values, test2(1).values))
hold on 
for i = 2:numel(test2)
    plot(test2(i).age*(-1)+1950,PERCENTRANK(test2(i).values, test2(i).values))
end
hold off



test5 = PERCENTRANK(test2(1).values, test2(1).values)

PERCENTRANK = @(YourArray, TheProbes) reshape( mean( bsxfun(@le, YourArray(:), TheProbes(:).') ) * 100, size(TheProbes) );