function result = runArcticAtlantic
% RUNARCTICATLANTIC Uses PaiCo to reconstruct the temperature signal for
%   the Arctic part of Atlantic.

mpath = fileparts(mfilename('fullpath'));
load([mpath filesep '..' filesep 'data' filesep 'arcticatlantic.mat']);

data.target.times = 00:2000;

options.resample.number = 000;
options.resample.parallelize = true;
options.filecache = false;
options.damping = [];
options.errorTolerance = 1e-8;
options.maxIters = 1e2;
options.heuristicstart = true;
options.regcov = 100;

tracker();
tic;
result = paico(data,options,@tracker);
toc;
save([mpath filesep 'results' filesep 'arcticatlantic.mat']);