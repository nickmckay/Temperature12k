function runPseudoProxy(datafile, methods)

load(datafile);
mpath = fileparts(mfilename('fullpath'));
if (nargin < 2)
    methods = {'PaiCo' 'MoM' 'OLS' 'PCReg' 'RegEM' 'LNA' 'BARCAST'};
    timeResult = cell(numel(dataCollection),1);
    pseudoResult = cell(numel(dataCollection),1);
    noiseResult = zeros(numel(dataCollection),2);    
else
    if (~iscell(methods)) 
        methods = {methods};
    end
    load([mpath filesep 'results' filesep 'pseudoexperiment.mat'],'timeResult','pseudoResult','noiseResult');
end

paicoopts.regcov = 100;
paicoopts.maxIters = 1e2;
paicoopts.errorTolerance = 1e-9;

regemopts.stagTolerance = 1e-6;
regemopts.maxIterations = 1e4;
regemopts.regpar = [];

lnaopts.preSamplerIterations = 400;
lnaopts.samplerIterations = 400;
lnaopts.useModes = true;
lnaopts.MHparams.phi = [0.01^2 40];

barcastopts.preSamplerIterations = 1000;
barcastopts.samplerIterations = 1000;
barcastopts.useModes = true;
barcastopts.useSpatialCache = true;
% First iterations have great belief in instrumental data
barcastopts.priors.tau2_I = [0.5 0.001];
barcastopts.sampleCompleteField = true;

pcregopts.nrPatterns = [];

disp(['Total data: ' num2str(numel(dataCollection))]);
ppm = ParforProgMon('Pseudoproxy ', numel(dataCollection), 1, 300, 80);
tic;
% Sort data so that estimating remaining time is more accurate
index = fastrandperm(1:numel(dataCollection));
dataCollection = dataCollection(index);
pseudoResult = pseudoResult(index);
timeResult = timeResult(index);

parfor i = 1:numel(dataCollection)
    data = dataCollection{i};
    result = pseudoResult{i};    
    times = timeResult{i};
    
    % Store target data
    result.target = data.target.data;

    if (any(strcmpi('paico',methods)))   
        % PaiCo
        tic;
        res = paico(data, paicoopts);
        times.PaiCo = toc;
        result.PaiCo = res.signal;
        noiseResult(i,:) = [data.target.noisestd res.noisestd];
    end

    if (any(strcmpi('mom',methods)))   
        % Method of Moments
        tic;
        res = mom(data);
        times.MOM = toc;
        result.MOM = res.signal;    
    end
    
    if (any(strcmpi('ols',methods)))   
        % Ordinary Least Squares
        tic;
        res = ols(data);
        times.OLS = toc;
        result.OLS = res.signal;
    end
    
    if (any(strcmpi('pcreg',methods)))       
        % Principal Component Regression
        tic;
        res = pcreg(data, pcregopts);
        times.PCReg = toc;
        result.PCReg = res.signal;
    end
    
    if (any(strcmpi('regem',methods)))   
        % RegEM with most significant eigenvector
        tic;
        res = regem(data, regemopts);
        times.RegEM = toc;
        result.RegEM = res.field;
    end
    
    if (any(strcmpi('lna',methods)))       
        % LNA
        tic;
        res = lna(data, lnaopts);
        times.LNA = toc;
        result.LNA = mean(res.signals,1);
    end

    if (any(strcmpi('barcast',methods)))   
        % BARCAST
        tic;
        res = barcast(data, barcastopts);
        times.BARCAST = toc;
        result.BARCAST = squeeze(mean(res.fields(:,:,500:end),3));     
    end
    
    pseudoResult{i} = result;
    timeResult{i} = times;
    
    ppm.increment();
end
[~,ind] = sort(index);
pseudoResult = pseudoResult(ind);
noiseResult = noiseResult(ind,:);
timeResult = timeResult(ind);

toc;
clear dataCollection;
ppm.delete()
mpath = fileparts(mfilename('fullpath'));
save([mpath filesep 'results' filesep 'pseudoexperiment.mat']);


