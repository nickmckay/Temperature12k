function calcAtlanticUncertainties

mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'results' filesep 'arcticatlantic.mat']);

samples = result.resample.signals;
nr = size(samples,1);

% Define time divisions
RWP = result.times >= 1 & result.times < 600;
DCAP = result.times >= 600 & result.times < 900;
MCA = result.times >= 900 & result.times < 1200;
LIA = result.times >= 1200 & result.times < 1900;
modern = result.times >= 1900;


% Mean of 1800-1900 below 0 and 1900-2000 above zero
prevyears = result.times >= 1800 & result.times < 1900;
curyears = result.times >= 1920 & result.times < 2000;
parts = [mean(samples(:,prevyears),2) < 0 ...
            (mean(samples(:,curyears),2) > 0)];
disp(['Pr(Mean 19st < 0)=' num2str(sum(parts(:,1))/nr)]);
disp(['Pr(Mean 20st > 0)=' num2str(sum(parts(:,2))/nr)]);

% Uncertainty for cooling trend
A = [ones(size(samples,2),1) result.times(:)];
alpha = A\(result.signal')*1e3;
disp(['Reconstruction 2kyr trend ' num2str(alpha(2,:))]);

alpha = A\samples'*1e3;
disp(['Ensemble 2kyr trend ' num2str(mean(alpha(2,:))) ' +- ' num2str(std(alpha(2,:)))]);

% Uncertainty in 2 deg variability
% Lower bound is the minimum of the higher value of each two consecutive years
% Upper bound the inverse
sr = sort(result.signal);
ss = sort(samples,2,'ascend');
dr = sr(:,round(size(sr,2)*.975)) - sr(:,round(size(sr,2)*.075));
dr = round(dr*10)/10;
ds = ss(:,round(size(ss,2)*.975)) - ss(:,round(size(ss,2)*.075));

disp(['Reconstruction: 95% of values within ' num2str(dr)]);
disp(['Ensemble: 95% of values within ' num2str(mean(ds)) '+-' num2str(std(ds))]);

% RWP heigher than any other
maxInstru = max(data.instrumental.data(data.instrumental.times <= 2000));
RWPtemp = sort(samples(:,RWP),2,'descend');
moderntemp = sort(samples(:,modern),2,'descend');
parts = [RWPtemp(:,5) max(moderntemp(:,1),maxInstru)];%ones(nr,1)*max(data.instrumental.data(:))];
disp(['RWP peaks highest p=' num2str(sum(parts(:,1) > parts(:,2))/nr)]);

% MCA equally high as modern
MCAtemp = sort(samples(:,MCA),2,'descend');
parts = [MCAtemp(:,5) max(moderntemp(:,1),maxInstru)];
disp(['MCA peaks highest p=' num2str(sum(parts(:,1) > parts(:,2))/nr)]);

% Comparisons to other reconstructions
pdata = xlsread('kaufman.xlsx');
ptimes = pdata(:,1);
pdata = pdata(:,2:end)';

pRWP = ptimes >= min(result.times(RWP)) & ptimes < max(result.times(RWP));
pMCA = ptimes >= min(result.times(MCA)) & ptimes < max(result.times(MCA));
RWPepisode = result.times >= 380 & result.times < 420;
pEpisode = ptimes >= 380 & ptimes < 420;
episode2 = result.times >= 1800 & result.times < 1880;
pEpisode2 = ptimes >= 1800 & ptimes < 1880;


rMeans = [mean(result.signal(MCA)) mean(result.signal(RWP)) ...
          mean(result.signal(RWPepisode)) mean(result.signal(episode2))];

pMeans = [mean(pdata(:,pRWP),2) mean(pdata(:,pMCA),2) ...
          mean(pdata(:, pEpisode), 2) mean(pdata(:, pEpisode2),2)];

sMeans = [mean(samples(:,MCA),2) mean(samples(:,RWP),2) ...
          mean(samples(:,RWPepisode),2) mean(samples(:,episode2),2)];

sampleHotter = [mean(sMeans-repmat(max(pMeans,[],1),size(sMeans,1),1),1);...
          std(sMeans-repmat(max(pMeans,[],1),size(sMeans,1),1),1)];
      
sampleCooler = [mean(sMeans-repmat(min(pMeans,[],1),size(sMeans,1),1),1);...
          std(sMeans-repmat(min(pMeans,[],1),size(sMeans,1),1),1)];

disp('Differences over episodes:');
disp(['RWP: Reconstruction ' num2str(rMeans(1) - max(pMeans(:,1))) ' hotter (' num2str(sampleHotter(1,1)) '+-' num2str(sampleHotter(2,1)) ')']);
disp(['MCA: Reconstruction ' num2str(rMeans(2) - max(pMeans(:,2))) ' hotter (' num2str(sampleHotter(1,2)) '+-' num2str(sampleHotter(2,2)) ')']);
disp(['RWP episode: Reconstruction ' num2str(rMeans(3) - max(pMeans(:,3))) ' hotter (' num2str(sampleHotter(1,3)) '+-' num2str(sampleHotter(2,3)) ')']);
disp(['Around 1840: Reconstruction ' num2str(rMeans(4) - max(pMeans(:,4))) ' hotter (' num2str(sampleCooler(1,4)) '+-' num2str(sampleCooler(2,4)) ')']);

mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'results' filesep 'hotcools.mat']);

hotcools(size(hotcools,1)+1,:) = [find(result.times == 1900) find(result.times==1940)];
for i = 1:size(hotcools,1)
    A = [ones(hotcools(i,2)-hotcools(i,1)+1,1) (hotcools(i,1):hotcools(i,2))'];
    alphaRes = A\(result.signal(hotcools(i,1):hotcools(i,2))')*10;
    alphaSam = A\samples(:,hotcools(i,1):hotcools(i,2))'*10;    
    disp(['Trend at ' num2str(result.times(hotcools(i,1))) '-' num2str(result.times(hotcools(i,2))) ...
         ' has slope ' num2str(alphaRes(2)) ' (' num2str(mean(alphaSam(2,:))) ' +- ' num2str(std(alphaSam(2,:))) ') C/10yrs']);
end

% Difference of averages of MCA and LIA
fLIA = result.times >= 1601 & result.times <= 1630;
fMCA = result.times >= 1071 & result.times <= 1100;

avgDiff = mean(samples(:,fMCA),2) - mean(samples(:,fLIA),2);
rDiff = mean(result.signal(fMCA)) - mean(result.signal(fLIA));
disp(['Frank et al. 2010 times for average difference of MCA and LIA ' num2str(rDiff) ' (' num2str(mean(avgDiff)) ' +- ' num2str(std(avgDiff)) ')']);

gLIA = result.times >= 1620 & result.times <= 1650;
avgDiff = mean(samples(:,fMCA),2) - mean(samples(:,gLIA),2);
rDiff = mean(result.signal(fMCA)) - mean(result.signal(gLIA));
disp(['Goosse et al. 2012 times for average difference of MCA and LIA ' num2str(rDiff) ' (' num2str(mean(avgDiff)) ' +- ' num2str(std(avgDiff)) ')']);

% RWP interval 50-150
RWPinterval = result.times >= 50 & result.times < 150;
RWPint = sort(samples(:,RWPinterval),2,'descend');
parts = [RWPint(:,5) moderntemp(:,1)];%ones(nr,1)*max(data.instrumental.data(:))];
disp(['RWP interval 50-150 peaks highest p=' num2str(sum(parts(:,1) > parts(:,2))/nr)]);

% RWP episode
RWPepi = sort(samples(:,RWPepisode),2,'descend');
parts = [RWPepi(:,5) moderntemp(:,1)];%ones(nr,1)*max(data.instrumental.data(:))];
disp(['RWP 390-420 peaks highest p=' num2str(sum(parts(:,1) > parts(:,2))/nr)]);

return;

% Average millenial scale cooling trend
slope = zeros(nr+1, size(samples,2)-1000 + 1);
A = [ones(1e3,1) (1:1e3)'];
ind = 0:(1e3-1);
for i = 1:size(slope,2)
    coeff = A\[result.signal(i+ind);samples(:,i+(ind))]';
    slope(:,i) = coeff(2,:);
end
slope = mean(slope,2)*1000;
disp(['Average cooling trend ' num2str(slope(1)) ' C/1000yrs (' num2str(mean(slope(2:end))) ' +- ' num2str(std(slope(2:end))) ')' ]);