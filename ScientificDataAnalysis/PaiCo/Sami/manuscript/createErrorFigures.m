function createErrorFigures(tracker)

if (~exist('tracker','var'))
    tracker = @(a,b,c)([]);
end

mpath = fileparts(mfilename('fullpath'));
load([mpath filesep 'manuscript_colors.mat']);
load([mpath filesep 'results' filesep 'pseudoexperiment.mat']);

% Convert data to more manageable format
errors = zeros(numel(nrproxy), numel(snr), runs*2, 129);

tracker('Error figures');
pos = 0;
for ni = 1:numel(nrproxy)
    for si = 1:numel(snr)
        stats = zeros(runs*2,129);        
        for nl = 1:2
            for ri = 1:runs
                pos = pos + 1;


                target = pseudoResult{pos}.target';
                paico = pseudoResult{pos}.PaiCo';
                
                St = pwelch(target);
                Sp = pwelch(paico);
                
                stats((nl-1)*runs + ri,:) = (Sp - St)/sum(St);
                errors(ni,si,(nl-1)*runs + ri,:) = (Sp - St)/sum(St);
                tracker('Error figures', pos, numel(pseudoResult));
            end     
        end
%         figure(1); clf; hold all;
%         mn = mean(stats,1);
%         sd = std(stats,[],1);
%         plot(mn);
%         plot(mn+sd,'k:');
%         plot(mn-sd,'k:');
%         plot([0 1],[0 0],'r:');
%         title([num2str(nrproxy(ni)) ' ' num2str(snr(si))]);
%         drawnow;
    end
end

mpath = fileparts(mfilename('fullpath'));
colors = [darkblue; brown; magenta; cyan];
figure(5); clf; hold all;
names = cell(numel(nrproxy)*2,1);
w = 1./linspace(0,.5,129);
w(1) = targetLength;
for ni = 1:numel(nrproxy)
    mask = snr < 1;
    plot(w,squeeze(mean(mean(errors(ni,mask,:,:),2),3)),':','color',colors(ni,:));
    plot(w,squeeze(mean(mean(errors(ni,~mask,:,:),2),3)),'-','color',colors(ni,:));
    names{(ni-1)*2+1} = ['M=' num2str(nrproxy(ni)) ', SNR < 1'];
    names{ni*2} = ['M=' num2str(nrproxy(ni)) ', SNR > 1'];
end
legend(names);
set(gca,'XScale','log','XDir','reverse','XTick',[2 5 10 20 50 100 200 500]);
axis tight;
print('-dpdf',[mpath filesep 'figures' filesep 'spectralerrors.pdf']);
