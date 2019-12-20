function [options, info] = defaults(opt)
% Defaults of PaiCo

info.method.name = 'Pairwise Comparisons';
info.method.description = '';
info.method.authors = 'Sami Hanhijärvi';
info.method.doi = '';
info.method.reference = '';

info.resample.number = ...
    ['Number of resampled results. Default = 0, i.e., no resampling.'];         
options.resample.number = setDefault(opt, {'resample' 'number'}, 0);

info.resample.parallelize = ...
    ['Boolean stating whether or not to parallelize resampling. Default = false. '];         
options.resample.parallelize = setDefault(opt, {'resample' 'parallelize'}, false);

info.fileCache = ...
    ['Boolean stating whether or not to use a flat file for storing massive matrices. ' ...
     'Cannot be used if resampling is parallelized. Default = false.'];
options.fileCache = setDefault(opt, {'fileCache'}, false);

info.regcov = ...
    ['Covariance of the regularization for the temperature signal. If matrix, has to be '...
    'the same length as target. A scalar is interpreted as diagonal covariance matrix ' ...
    'with given scalar on each diagonal element. Default = 10. Required.'];
options.regcov = setDefault(opt, {'regcov'}, 10);

info.heuristicStart = ...
    ['Boolean stating whether or not to being the maximum likelihood search using a ' ...
    'heuristic initial value. The initial value is calculated by summing the rows of ' ...
    'linear comparisons matrix. Default = true.'];
options.heuristicStart = setDefault(opt, {'heuristicStart'}, true);

info.errorTolerance = ...
    ['Scalar value depicting the maximum relative error tolerance in the ML iteration.' ...
    'Increasing this value will halt the iteration quicker, but may produce suboptimal results.' ...
    'Default = 1e-9.'];
options.errorTolerance = setDefault(opt, {'errorTolerance'}, 1e-9);

info.maxIters = ...
    ['Maximum number of iterations carried out in ML optimization. Default = 1e4.'];
options.maxIters = setDefault(opt, {'maxIters'}, 1e4);

info.randomStarts = ...
    ['The number of times the initial signal is randomized and optimization is carried out again. ' ...
    'Has no effect if heuristicStart is set to true. Default = 1.'];
options.randomStarts = setDefault(opt, {'randomStarts'}, 1);


info.damping = ...
    ['Scalar multipier used for Levenberg-Marquardt damping. Empty equals adaptive. ' ...
    'Default = [].'];
options.damping = setDefault(opt, {'damping'}, []);

if (options.resample.parallelize)
    parainfo = gcp;

    if (~exist('parpool'))
        error('RECON:PAICO','Parallelization of resampling requested but Parallel Computing Toolbox not found.');
    end
    if (options.fileCache)
        error('RECON:PAICO','Both file cache and resampling parallelization cannot be active at the same time.');
    end

    if (parainfo.NumWorkers <= 1)
        warning('RECON:PAICO','parpool is 1 or less. Please initialize parpool to at least 2 to parallelize.');
    end
end

function val = setDefault(opt, field, default)

val = default;
for i = 1:numel(field),
    if (isfield(opt, field{i}))
        opt = opt.(field{i});        
    else
        return;
    end
end
val = opt;