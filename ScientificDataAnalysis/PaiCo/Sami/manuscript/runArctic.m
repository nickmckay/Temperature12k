% function result = runArctic
% % RUNARCTICATLANTIC Uses PaiCo to reconstruct the temperature signal for
% %   the Arctic part of Atlantic.

% mpath = fileparts(mfilename('fullpath'));
% load([mpath filesep '..' filesep 'data' filesep 'arcticatlantic.mat']);
%%
load OrigRecon

arc.target.times = 00:2000;
arc.instrumental.noisestd = .2224;
arc.instrumental.data = OrigRecon.Arc.Temp(1851:end)';
arc.instrumental.times = OrigRecon.Arc.Year(1851:end)';

options.resample.number = 100;
options.resample.parallelize = true;        
options.filecache = false;  
options.damping = [];
options.errorTolerance = 1e-8;
options.maxIters = 1e2;
options.heuristicstart = true;
options.regcov = 1000;
delete(gcp)
matlabpool local 4
addpath /home/nick/Dropbox/ML_Scripts/Sami/common
pctRunOnAll javaaddpath /home/nick/Dropbox/ML_Scripts/Sami/common
tracker();
tic;
result = paico(arc,options,@tracker);
toc;
%save([mpath filesep 'results' filesep 'arctic.mat']);

%%
%estimate 95% cl
s_ens=sort(result.resample.signals);
cl95=s_ens(ceil(.95*size(s_ens,1)),:);
cl05=s_ens(ceil(.05*size(s_ens,1)),:);

load OrigRecon

close all
plot(OrigRecon.Arc.Year,OrigRecon.Arc.Temp,'k')
hold on
plot(OrigRecon.Arc.Year,OrigRecon.Arc.Min,'k--','color',[.7 .7 .7 ])
plot(OrigRecon.Arc.Year,OrigRecon.Arc.Max,'k--','color',[.7 .7 .7 ])

plot(result.times,result.signal,'r')
plot(result.times,cl95,'b')
plot(result.times,cl05,'b')


