function createPseudoFigures



mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'results' filesep 'pseudoexperiment.mat']);
mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'manuscript_colors.mat']);

alpha = 0.05;


% Convert data to more manageable format
names = fieldnames(pseudoResult{1});
targetInd = find(strcmpi(names, 'target'));
methodInd = setdiff(1:numel(names),targetInd);

% results = nan(numel(nrproxy), numel(snr), 2, runs, numel(methodInd), targetLength);
stats = nan(numel(nrproxy), numel(snr), 2, numel(methodInd), 4);
meanranks = nan(numel(nrproxy), numel(snr), 2, numel(methodInd), 4);
meantimes = nan(numel(nrproxy), numel(snr), 2, numel(methodInd));
worseThanBest = false(numel(nrproxy), numel(snr), 2, numel(methodInd), 4);
noises = zeros(numel(nrproxy), numel(snr), 2, runs, 2);

pos = 0;
for ni = 1:numel(nrproxy)
    for si = 1:numel(snr)
%         coverage = zeros(1,numel(coverAlpha));
        for nl = 1:2
            runstats = zeros(numel(methodInd),4, runs);
            times = zeros(numel(methodInd),1);
            for ri = 1:runs
                pos = pos + 1;
                noises(ni,si,nl,ri,:) = noiseResult(pos,:);

                target = pseudoResult{pos}.target;

                % Calculate statistics                
                for mi = 1:numel(methodInd)
                    signal = pseudoResult{pos}.(names{methodInd(mi)});
                    times(mi) = times(mi) + timeResult{pos}.(names{methodInd(mi)});
                    mn = mean(target((numel(target) - calibLength + 1):end));
                    runstats(mi,:,ri) = [sqrt(mean((signal - target).^2))/std(target) ... % RMSE
                        corr(signal(:), target(:)) ... % Correlation
                        (1-sum((target-signal).^2)/sum((target-mn).^2)) ... % Reduction in Error
                        (1-sum((target-signal).^2)/sum((target-mean(target)).^2))]; % Coefficient of Efficiency                    
                end
            end                 
                        
            meantimes(ni, si, nl, :) = times/runs;
            stats(ni, si, nl, :, :) = squeeze(mean(runstats,3));            
            % For RMSE, smaller is better, and for others, larger is
            % better. Flip RMSE to get comparable rankings.
            runstats(:,1,:) = -runstats(:,1,:);             
            [~, inds] = max(squeeze(mean(runstats,3)),[],1);
            for ii = 1:4
                % Limit values to the worst that the best method had. Very
                % poorly performed methods have a large variance and the
                % following t-test can not reject the null hypothesis
                % because of that. Limiting the values yields figures that 
                % better reflect the performance.
                s = squeeze(runstats(:,ii,:));
                limit = min(s(inds(ii),:));
                s(s < limit) = limit;
                
                worseThanBest(ni, si, nl, :, ii) = ttest2(s', repmat(squeeze(runstats(inds(ii), ii, :)),1, size(runstats,1)), alpha, 'left');
            end
            [~, inds] = sort(runstats,1,'descend');
            ranks = zeros(size(inds));
            for i = 1:size(inds,2)
                for ri = 1:runs
                    ranks(sub2ind(size(ranks), ...
                        squeeze(inds(:,i,ri)), ...
                        ones(size(inds,1),1)*i, ...
                        ones(size(inds,1),1)*ri)) = 1:numel(methodInd);
                end
            end
            meanranks(ni, si, nl, :, :) = squeeze(mean(ranks,3));
        end
%         coverages(ni,si,:) = coverage/(runs*2);
    end
end

% for si = 1:4
%     for ni = 1:numel(nrproxy)
%         for nl = 1:2
%             figure(4); clf; hold all;
%             plot(squeeze(stats(ni, :, nl, :, si)));
%             if (si == 1)
%                 plot(1./snr/sqrt(nrproxy(ni)),'r:');
%             end
%             legend(names(methodInd));
%             title(['n=' num2str(nrproxy(ni)) ' linear=' num2str(nl-1) ' stat=' num2str(si)]);
%             set(gca,'XTick',1:numel(snr));            
%             set(gca,'XTickLabel',num2str(snr'));
%             axis([1 numel(snr) 0 1]);
%             
%             pause;
%         end
%     end
% end

mpath = fileparts(mfilename('fullpath'));

methodColors = [darkblue; green; cyan; brown; magenta; lightgray; lightblue];

% Execution times
figure(1); clf; 
bar((squeeze(mean(mean(mean(meantimes, 1), 2), 3))));
set(gca,'XTickLabel', names(methodInd));
squarepage([8,8]);
print('-dpdf',[mpath filesep 'figures' filesep 'meantimes.pdf']);

% PaiCo's ability to recover noise variance
figure(2); clf; hold all;
plot(noiseResult(:,1), noiseResult(:,2), '.');
mx = max(noiseResult(:));
axis([0 mx 0 mx]);
plot([0 mx], [0 mx]);
text(0.1,mx*.8,['RMSE=' num2str(sqrt(mean((noiseResult(:,1) - noiseResult(:,2)).^2)))]);
text(0.1,mx*.7,['R^2=' num2str(corr(noiseResult(:,1), noiseResult(:,2))^2)]);
squarepage([8,8]);
print('-dpdf', [mpath filesep 'figures' filesep 'noiserecovery.pdf']);

% Plot ranks
figure(3); clf;
meanranks = permute(meanranks,[4 5 3 2 1]);
meanranks = squeeze(mean(meanranks(:,:,:,:), 4));
subplot(2,1,1); hold all;
bar(squeeze(meanranks(:,:,1))');
subplot(2,1,2); hold all;
bar(squeeze(meanranks(:,:,2))');
legend(names(methodInd),'location','eastoutside');
squarepage([8,8]);
print('-dpdf', [mpath filesep 'figures' filesep 'pseudoranks.pdf']);
            
% Plot different statistics for different parameters
figure(4); clf; hold all;
ni = 1; nl = 1; si = 1;
for mi = 1:size(stats,4)
    plot(squeeze(stats(ni, :, nl, mi, si)),'color',methodColors(mi,:));
end
plot(1./snr/sqrt(nrproxy(ni)),'r:');
legend(names(methodInd));
title(['n=' num2str(nrproxy(ni)) ' linear=' num2str(nl-1) ' stat=' num2str(si)]);
set(gca,'XTick',1:numel(snr));            
set(gca,'XTickLabel',num2str(snr'));
axis([1 numel(snr) 0 1]);
squarepage([8,8]);
print('-dpdf', [mpath filesep 'figures' filesep 'pseudo_10_nonlin_rmse.pdf']);

figure(5); clf; hold all;
ni = 3; nl = 2; si = 1;
colormap([methodColors; 1 0 0]);
for mi = 1:size(stats,4)
    plot(squeeze(stats(ni, :, nl, mi, si)),'color',methodColors(mi,:));
end
plot(1./snr/sqrt(nrproxy(ni)),'r:');
legend(names(methodInd));
title(['n=' num2str(nrproxy(ni)) ' linear=' num2str(nl-1) ' stat=' num2str(si)]);
set(gca,'XTick',1:numel(snr));            
set(gca,'XTickLabel',num2str(snr'));
axis([1 numel(snr) 0 1]);
squarepage([8,8]);
print('-dpdf', [mpath filesep 'figures' filesep 'pseudo_50_lin_rmse.pdf']);

% Plot worstThanBest statistics
figure(6); clf; hold all;
line([1 5], [1 1]);
line([1 5], [2 2]);
radius = 0.02;
% Only include PaiCo, MoM, LNA and BARCAST as OLS, PCReg and RegEM were
% always worse than the others
includeInd = [1 2 6 7]; 
squarepage([8,8]);
for si = 1:4
    line([si si], [1 3]);
    for nl = 1:2
        yy = linspace(0, 1, numel(snr)+2);
        xx = linspace(0, 1, numel(nrproxy)+2);
        yy = nl + repmat(yy(2:(end-1))',1,numel(nrproxy));
        xx = si + repmat(xx(2:(end-1)), numel(snr), 1);
        plot(xx,yy,'r.');
        if (si == 1)
            text(ones(numel(snr),1),yy(:,1), num2str(snr'),'HorizontalAlignment','right');
        end
        if (nl == 1)
            text(xx(1,:)', ones(numel(nrproxy),1), num2str(nrproxy'));
        end
        for mi = 1:numel(includeInd)
            mask = squeeze(~worseThanBest(:,:,nl,includeInd(mi),si))';
            plot(xx(mask) + cos(mi/(numel(includeInd))*2*pi+pi/4)*radius, ...
                yy(mask) + sin(mi/(numel(includeInd))*2*pi+pi/4)*radius, ...
                '.','color',methodColors(includeInd(mi),:),'markerSize',16);
        end
    end
end
axis off;

print('-dpdf', [mpath filesep 'figures' filesep 'pseudo_best_methods.pdf']);
