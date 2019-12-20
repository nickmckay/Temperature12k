function [params, signal] = initialValues(data, model, options)
%% Initializes the parameter values for the MCMC procedure 

if options.useModes
    % Priors for proxies
    for i = 1:numel(data.proxy)
        params.betaL{i} = ones(size(data.proxy{i}.data,1),1)*options.priors.betaL(i,1);
        params.muL{i} = ones(size(data.proxy{i}.data,1),1)*options.priors.muL(i,1);
    end
    params.beta = options.priors.beta(:,1);
    params.sigma2L = options.priors.sigma2L(:,2)./(options.priors.sigma2L(:,1)+1);
    params.phi2L = mean(options.priors.phi2L,2);
    params.phi1L = zeros(size(params.phi2L));    
else
    for i = 1:numel(data.proxy)
        params.betaL{i} = rand(size(data.proxy{i}.data,1),1)*sqrt(options.priors.betaL(i,2)) + options.priors.betaL(i,1);
        params.muL{i} = rand(size(data.proxy{i}.data,1),1)*sqrt(options.priors.muL(i,2)) + options.priors.muL(i,1);
    end
    
    params.beta = rand(size(options.priors.beta,1),1).*sqrt(options.priors.beta(:,2)) + options.priors.beta(:,1);
    params.phi2L = rand(size(options.priors.phi2L,1),1).*(diff(options.priors.phi2L,1,2))+options.priors.phi2L(:,1);
    
    params.phi1L = rand(size(options.priors.phi2L,1),1).*(diff(options.priors.phi2L,1,2))+options.priors.phi2L(:,1);
    params.phi1L = params.phi1L.*(options.priors.phi2L(:,2)-params.phi2L);
          
    %Initial Value for sigma2's: draw from inverse gamma prior is likely a bad
    %ideas, due to the very large possible values. So truncate to below some
    %value:    
    params.sigma2L = zeros(size(options.priors.sigma2L,1),1);
    for i = 1:size(options.priors.sigma2L,1)        
        while (1)
            t = 1/gamrnd(options.priors.sigma2L(i,1),1/options.priors.sigma2L(i,2));
            if (t < 5)
                params.sigma2L(i) = t;
                break;
            end
        end    
    end
end

if options.useModes
    signal = nanmean(data.instrumental.data(:))*ones(1,numel(model.timeline));    
else
    signal = nanmean(data.instrumental.data(:)) + randn(size(model.timeline))*nanstd(data.instrumental.data(:));
end