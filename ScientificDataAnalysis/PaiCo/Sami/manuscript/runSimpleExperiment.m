function runSimpleExperiment(redo)

runs = 40;

mpath = fileparts(mfilename('fullpath'));
storefile = [mpath filesep 'results' filesep 'simpleexperiment_new.mat'];
if (~exist('redo','var') || ~exist(storefile,'file'))
    redo = 1:5;
else
    % Redo portion of experiments
    load(storefile,'results','jobs');
end

defaultJob.nrproxy = 20;
defaultJob.noisevar = 1;
defaultJob.param = [];
defaultJob.startHeight = 10;
defaultJob.transfer = @(x)(x);

if (ismember(1,redo))
    job = defaultJob;
    job.noisevar = .5;
    jobs(1) = job;
end

if (ismember(2,redo))
    job = defaultJob;
    p = [0.1 0.25 0.5 0.75 1 1.5 2 3];
    job.noisevar = repmat(p',1,job.nrproxy);
    job.param = {p};
    jobs(2) = job;
end

if (ismember(3,redo))
    job = defaultJob;
    job.nrproxy = [5 10 20 30 40 60 80 100]';
    job.param = job.nrproxy;
    jobs(3) = job;
end

if (ismember(4,redo))
    job = defaultJob;
    c = linspace(-1,1,job.nrproxy);
    p = 0:0.1:1;
    job.param = p;
    job.noisevar = exp(-repmat(p(:),1,numel(c)).*repmat(c(:)',numel(p),1));
    jobs(4) = job;
end

if (ismember(5,redo))
    job = defaultJob;      
    p = [0.01 0.02:0.02:0.2];
    job.param = p;
    job.transfer = cell(numel(p),1);
    for i = 1:numel(p)
        job.transfer{i} = @(x)(atan(p(i)*x));
    end
    jobs(5) = job;
end

for i = redo(:)'
    results{i} = analyze(jobs(i), runs, i);
    disp(['Job ' num2str(i) ' done']);
    save(storefile,'results','jobs','runs');
end


function jobResult = analyze(batchjob, runs, jobind)

paicoopts.regcov = 100;
paicoopts.errorTolerance = 1e-8;
paicoopts.maxIters = 1e2;

regemopts.stagTolerance = 1e-6;
regemopts.maxIterations = 1e4;
regemopts.regpar = [];

lnaopts.preSamplerIterations = 400;
lnaopts.samplerIterations = 300;
lnaopts.useModes = true;
lnaopts.MHparams.phi = [0.01^2 20];

barcastopts.preSamplerIterations = 1000;
barcastopts.samplerIterations = 1000;
barcastopts.useModes = true;
barcastopts.useSpatialCache = true;
% First iterations have great belief in instrumental data
barcastopts.priors.tau2_I = [0.5 0.001];
barcastopts.sampleCompleteField = true;

pcregopts.nrPatterns = 1;

tic;

fields = fieldnames(batchjob);
nn = zeros(numel(fields),1);
for i = 1:numel(fields)
    nn(i) = size(batchjob.(fields{i}),1);
end

[count,ind] = max(nn);
fieldind = ind(1);
count = count(1);

ppm = ParforProgMon(['Simple: job ' num2str(jobind)], runs*count);

jobResult  = cell(count, runs);
for ti = 1:count
    job = batchjob;
    t = job.(fields{fieldind});
    job.(fields{fieldind}) = t(ti,:);
    
    if (iscell(job.transfer))
        job.transfer = job.transfer{1};
    end
    
    if (numel(job.noisevar) == 1)
        job.noisevar = ones(job.nrproxy,1)*job.noisevar;
    end
    
    parfor i = 1:runs
        data = [];
        result = [];

        data.instrumental.times = 401:500;
        t = [job.startHeight*ones(1,100) linspace(job.startHeight,0,100) zeros(1,100) sin((1:200)/100*2*pi)];
        result.target = t;
        data.instrumental.data = t(data.instrumental.times);
        data.target.times = 1:500;
        data.target.data = t;
        data.proxy = cell(job.nrproxy,1);
        for j = 1:job.nrproxy
            data.proxy{j}.data = zscore(job.transfer(t + randn(1, numel(t))*job.noisevar(j)));
            data.proxy{j}.times = 1:500;
            data.proxy{j}.locations = [1 1];
        end
        data.instrumental.locations = [1 1];

        % PaiCo
        res = paico(data, paicoopts);
        result.PaiCo = res.signal;

        % Method of moments
        res = mom(data);
        result.MOM = res.signal;

        % Ordinary Least Squares
        res = ols(data);
        result.OLS = res.signal;    

        % Principal Component Regression
        res = pcreg(data, pcregopts);
        result.PCReg = res.signal;

        % RegEM with most significant eigenvector
        res = regem(data, regemopts);
        result.RegEM = res.field;  

        % LNA
        res = lna(data, lnaopts);
        result.LNA = mean(res.signals,1);

        % BARCAST  
        res = barcast(data, barcastopts);
        result.BARCAST = squeeze(mean(res.fields(:,:,500:end),3));

        jobResult{ti,i} = result;

        ppm.increment();
    end
end

toc;
ppm.delete();
