function [options, info] = defaults(opt, data)
%% Assign priors and MCMC default values for the given data

info.method.name = 'BARCAST';
info.method.description = ...
    'Bayesian Hieracrchical model for climate reconstructions.';
info.method.authors = 'Martin P. Tingley, Sami Hanhijärvi';
info.method.reference = ...
    ['Tingley, M, Huybers, P. A Bayesian Algorithm for Reconstructing ' ...
    'Climate Anomalies in Space and Time. Part I: Development and ' ...
    'Applications to Paleoclimate Reconstruction Problems. Journal of ' ...
    'Climate, 2010.'];
info.method.doi = '10.1175/2009JCLI3015.1';

info.priors.alpha = ...
    ['Prior distribution parameters for the autocorrelation of temperature. ' ...
    'Defines the range for uniform distribution.'];
options.priors.alpha = setDefault(opt, {'priors' 'alpha'}, [0,1]); %[min, max) of uniform.

info.priors.mu = ...
    ['Prior distribution paramters for the mean of temperature field. ' ...
    'Defines the mean and variance of gaussian distribution.'];
if (~isempty(data))
    options.priors.mu = setDefault(opt, {'priors' 'mu'},[nanmean(data.instrumental.data(:)), 5^2]);
else
    options.priors.mu = setDefault(opt, {'priors' 'mu'},[0, 5^2]);
end


info.priors.sigma2 = ...
    ['Prior distribution patameters for the partial sill of the spatial ' ...
    'covariance matrix of the innovations that drive the AR(1) process. ' ...
    'Defines the shift and scale, resp., of inverse gamma distribution.'];
options.priors.sigma2= setDefault(opt, {'priors' 'sigma2'}, [0.5,0.5]);
% This is equivalent to 1 observations with an average square deviation of 2
% Very general, broad prior. 
% note, calling the pars a and b, and the pars of the scaled inverse chi^2
% nu and s^2, we have nu=2a and s^2=b/a. These give, in the bayesian sense,
% the number of prior samples and the mean squared deviation.


info.priors.phi = ...
    ['Prior distribution paramaters for the inverse range of this spatial ' ...
    'covariance matrix. Defines the log-mean and log-variance of the ' ...
    'lognormal distribution.'];    
options.priors.phi=setDefault(opt, {'priors' 'phi'}, [-4.65, 1.2]);
%RECALL we are using a log transformation in the Metropolis sampler for
%this step, so the prior is log normal. There are nine data points evenly
%space by about 111km. A largely uninformative prior would specify the
%range to be somewhere between 10 and 1000 km, so the log of the 
%inverse range should be between -7 and -2.3

info.priors.tau2_I = ...
    ['Prior distribution parameters for the error variance of instrumental ' ...
    'observations. Defines the shift and scale, resp., of inverse gamma distribution.'];
options.priors.tau2_I=setDefault(opt, {'priors' 'tau2_I'}, [0.5,0.5]);
%This is equivalent to 1 observations (nu=2*a) with an average square
%deviation of 1 (s^2=a/b)

info.priors.tau2_P = ...
    ['Prior distribution parameters for the error variance of proxy ' ...
    'observations. Defines the shift and scale, resp., of inverse gamma distribution.'];
if (~isempty(data))
    options.priors.tau2_P = setDefault(opt, {'priors' 'tau2_P'}, ones(numel(data.proxy),1)*[0.5,0.5]);
else
    options.priors.tau2_P = setDefault(opt, {'priors' 'tau2_P'}, [0.5, 0.5]);
end
%This is equivalent to 1 observations (nu=2*alpha) with an average square
%deviation of 2


info.priors.Beta_1 = ...
    ['Prior distribution parameters for the scaling of proxy observations. '...
    'Defines the mean and variance of normal distribution.'];
%This one is normal as well. IF THE PROXIES HAVE BEEN PRE-PROCSSES TO HAVE
%MEAN=0, STD=1, THEN Hypothetically, the scaling should be
%(1-tau_P^2)(1-alpha^2)/sigma^2)^(+1/2). So set the mean to the modes of
%these priors, and then include a decent sized variance.
if ~isempty(data)
    options.priors.Beta_1= setDefault(opt, {'priors' 'Beta_1'}, [((1-options.priors.tau2_P(:,2)./(options.priors.tau2_P(:,1)+1))*(1-mean(options.priors.alpha)^2)/(options.priors.sigma2(2)/(options.priors.sigma2(1)+1))).^(1/2), ones(numel(data.proxy),1)*8^2]);
else
    options.priors.Beta_1 = setDefault(opt, {'priors' 'Beta_1'}, [1, 8^2]);
end

info.priors.Beta_0 = ...
    ['Prior distribution parameters for the shift of proxy observations. ' ...
    'Defines the mean and variance of normal distribution. Should be equal ' ...
    'to the prior of mu.'];
%Prior for Beta_0. SET EQUAL TO THE PRIOR FOR MU
if ~isempty(data)
    options.priors.Beta_0= setDefault(opt, {'priors' 'Beta_0'}, [-options.priors.Beta_1(:,1)*options.priors.mu(1) , ones(numel(data.proxy),1)*8^2]);
else
    options.priors.Beta_0=setDefault(opt, {'priors' 'Beta_0'}, [options.priors.mu(1), 8^2]);
end

info.priors.T_0 = ...
    ['Prior distribution parameters for the temperature field. ' ...
    'Defines the mean and variance of normal distribution.'];
if ~isempty(data)
    options.priors.T_0 = setDefault(opt, {'priors' 'T_0'}, [nanmean(data.instrumental.data(:)), 4*nanvar(data.instrumental.data(:))]);
else
    options.priors.T_0 = setDefault(opt, {'priors' 'T_0'}, [0, 8^2]);
end

%with mean given by mean of all inst and standard deviation 2
%times the std of all Inst.
%So the prior parameters (mean, var) are: 
%Note that this assumes a constant mean and diagonal cov mat. 

info.MHpars.log_phi = ...
    ['Parameters for the Metropolis-Hastings algorithm when sampling log phi.' ...
    'First is the variance of proposal distribution, and second is the number ' ...
    'of steps in the MH-sampling.'];
% Also set the MCMC jumping parameters
%The phi step.
% The jumping distribution of log(phi) is normal with mean zero; this sole paramter is
% the VARIANCE. We expect the posterior variance to be far smaller than the
% prior variance, so the jumping variance is set low. this can be adjusted
% as needed. 
%ALSO Specify the number of iterations of the Metropolis step for each iteration
%of the Gibbs: the Metroplis sampler needs to (come close to) converging each
%time. Until the algorithm settles down, this will take a little while.
%NOTE that this step is not time consuming. 
%First par is variance, second is number. 
options.MHpars.log_phi=setDefault(opt, {'MHpars' 'log_phi'}, [.04^2, 100]);

info.samplerIterations = ...
    ['Number of consecutive samples from the bayesian model.'];
options.samplerIterations=setDefault(opt, {'samplerIterations'}, 2000);

info.preSamplerIterations = ...
    ['Number of times to update only the temperature array before ' ...
    'beginning to update the other parameters.'];
options.preSamplerIterations=setDefault(opt, {'preSamplerIterations'}, 500); 

info.useModes = ...
    ['Use modes of priors as inital values for MCMC sampling.'];
options.useModes = setDefault(opt, {'useModes'}, false);

info.useSpatialCache = ...
    'Store spatial covariance matrices for speedy processing. Uses a lot of memory';
options.useSpatialCache = setDefault(opt, {'useSpatialCache'}, false);

info.sampleCompleteField = [...
    'Whether to sample the complete temperature field at once or sample ' ...
    'time separately. Increases convergence but requires massive amounts ' ...
    'memory. If BARCAST exits with "Out of memory" or "Maximum variable ' ...
    'size allowed by the program is exceeded.", this option should be set to false.' ...
    'Default = false.'];
options.sampleCompleteField = setDefault(opt, {'sampleCompleteField'}, false);

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