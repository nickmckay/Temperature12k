function createUncertFigures

mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'results' filesep 'uncertaintyExperiment_ols.mat']);
mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'manuscript_colors.mat']);

% Convert data to more manageable format
uncertainties = nan(numel(nrproxy), numel(snr), runs, targetLength);

pos = 0;
figure(1); clf; hold all;
colors = jet(numel(nrproxy)*numel(snr));

for ni = 1:numel(nrproxy)
    for si = 1:numel(snr)
        uncert = zeros(runs,targetLength);
        for ri = 1:runs
            pos = pos + 1;                
            uncert(ri,:) = sort(uncertResult{pos}.uncert);                 
            uncertainties(ni,si,ri,:) = uncert(ri,:);
%             plot(sort(uncertResult{pos}.uncert), (1:targetLength)/targetLength);
        end                                   
        mn = mean(uncert,1);
        sd = std(uncert,[],1);
        color = colors((ni-1)*numel(snr)+ si,:);
        plot(mn,(1:targetLength)/targetLength,'-','color',color);
        plot(mn-sd,(1:targetLength)/targetLength,':','color',color);
        plot(mn+sd,(1:targetLength)/targetLength,':','color',color);        
    end
end
plot([0 1],[0 1],'k-');
axis([0 1 0 1]);

figure(2); clf; hold all;
up = reshape(permute(uncertainties,[4 1 2 3]),[targetLength, numel(nrproxy)*numel(snr)*runs]);
mn = mean(up,2);
sd = std(up,[],2);

    
x = linspace(0,1,targetLength);
plot(x, (1:targetLength)'/targetLength - mn ,'k-');
plot(x, (1:targetLength)'/targetLength - (mn-sd),'k-');
plot(x,  (1:targetLength)'/targetLength - (mn+sd),'k-');

axis tight;

mpath = fileparts(mfilename('fullpath'));
% print('-dpdf',[mpath filesep 'figures' filesep 'uncertainty.pdf']);