function [options, info] = defaults(data, opt)
%% Assign priors and MCMC default values for the given data

if (~exist('data','var'))
    data= [];
end
if (~exist('opt','var'))
    opt = [];
end

info.method.name = 'LNA';
info.method.description = ...
    ['Bayesian Hieracrchical model for climate reconstructions. Does not yet' ...
    ' support boreholes.'];
info.method.authors = 'Li, Nychka and Ammann.';
info.method.reference = ['Li B, Nychka DW, Ammann CM (2010) The value of ' ...
    'multiproxy reconstruction of past climate. Journal of the American ' ...
    'Statistical Association  105(491):883–895.'];
info.method.doi = '10.1198/jasa.2010.ap09379';

nrProxy = 1;
if (~isempty(data))
    nrProxy = numel(data.proxy);
end

% Priors for proxies
info.priors.betaL = ...
    ['Prior distribution parameters for the scaling of proxy observations. '...
    'Defines the mean and variance of normal distribution.'];
%options.priors.betaL = setDefault(opt, {'priors' 'betaL'}, ones(nrProxy+1,1)*[1, 8^2]);
%NM change:
options.priors.betaL = setDefault(opt, {'priors' 'betaL'}, ones(nrProxy+1,1)*[0, 1]);
info.priors.muL = ...
    ['Prior distribution parameters for the shift of proxy observations. '...
    'Defines the mean and variance of normal distribution.'];
options.priors.muL = setDefault(opt, {'priors' 'muL'}, ones(nrProxy+1,1)*[0, 8^2]);

info.priors.beta = ...
    ['Prior distribution parameters for the forcing coefficients. ' ...
    'Defines the mean and variance of normal distribution.'];
if (~isempty(data) && isfield(data,'forcings'))
    options.priors.beta = setDefault(opt, {'priors' 'beta'}, ones(numel(data.forcings)+1,1)*[1, 8^2]);
else
    options.priors.beta = setDefault(opt, {'priors' 'beta'}, [1, 8^2]);
end

info.priors.sigma2L = ...
    ['Prior distribution parameters for the autoregressive noise variance. ' ...
    'Defines the shape and scale parameter of the inverse gamma distribution.'];
options.priors.sigma2L = setDefault(opt, {'priors' 'sigma2L'}, ones(nrProxy+1,1)*[0.5,0.5]);

info.priors.phi2L = ...
    ['Prior distribution parameters for the autoregressive noise lag-2 correlation. ' ...
    'Defines the range of uniform distribution.'];
%options.priors.phi2L = setDefault(opt, {'priors' 'phi2L'}, ones(nrProxy+1,1)*[-1,1]);
%NM change
options.priors.phi2L = setDefault(opt, {'priors' 'phi2L'}, ones(nrProxy+1,1)*[0,0]);

info.MGparams.phi = ...
    ['Parameters for the Metropolis-Hastings algorithm when sampling phis.' ...
    'First is the variance of proposal distribution, and second is the number ' ...
    'of steps in the MH-sampling.'];
options.MHparams.phi = setDefault(opt, {'MHparams' 'phi'}, [.01^2, 40]);

info.samplerIterations = ...
    ['Number of consecutive samples from the bayesian model.'];
options.samplerIterations = setDefault(opt, {'samplerIterations'}, 2000);

info.preSamplerIterations = ...
    ['Number of times to update only the temperature array before ' ...
    'beginning to update the other parameters.'];
options.preSamplerIterations = setDefault(opt, {'preSamplerIterations'}, 500);

info.useModes = ...
    ['Use modes of priors as inital values for MCMC sampling.'];
options.useModes = setDefault(opt, {'useModes'}, false);

if (~isempty(data))
    % Check data if given
    for i = 1:numel(data.proxy)
        if (isfield(data.proxy{i},'transform') && (numel(data.target.times) ~= size(data.proxy{i}.transform,1) || ...
                size(data.proxy{i}.transform,2) ~= size(data.proxy{i}.data,2)))
            error('RECON:LNA',['Transformation for proxy ' num2str(i) ' is of wrong size.']);
        end
    end

    if (isfield(data,'forcings'))
        for i = 1: numel(data.forcings)
            if ~isequal(data.forcings.times, data.target.times)
                error('RECON:LNA', ['Forcings need to be in the same times than the target. ' ...
                    'Transformations are not yet supported.']);
            end
        end
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