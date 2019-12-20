function [currentParams, currentField] = initialValues(data, model, options)
%% Initializes the parameter values for the MCMC procedure 
% Draws from (at times) truncated priors. 

if options.useModes
    %Initial Value of alpha: (uniform, so take mean):
    currentParams.alpha=mean(options.priors.alpha);

    %Initial Value for mu: mode of the normal prior.
    currentParams.mu=options.priors.mu(1);

    %Initial Value for sigma2: mode of inverse gamma prior:
    currentParams.sigma2=options.priors.sigma2(2)/(options.priors.sigma2(1)+1);

    %Initial Value for phi: mode of the log-normal distrubution:
    currentParams.phi=exp(options.priors.phi(1)-options.priors.phi(2));

    %Initial Value for tau2_I: : mode of inverse gamma prior:
    currentParams.tau2_I=options.priors.tau2_I(2)/(options.priors.tau2_I(1)+1);

    %Initial Value for tau2_P: : mode of inverse gamma prior:
    currentParams.tau2_P = options.priors.tau2_P(:,2)./(options.priors.tau2_P(:,1)+1);
    
    %Initial value for Beta_1: mode of the normal prior:
    currentParams.Beta_1 = options.priors.Beta_1(:,1);

    %Initial value for Beta_0: mode of the normal prior:
    currentParams.Beta_0 = options.priors.Beta_0(:,1);        

else

    %Initial Value of alpha: Draw from the uniform prior:
    %currentParams.alpha=rand(1)*diff(options.priors.alpha)+options.priors.alpha(1);
    %TRUCNATE A Bit to set initial value near 0.5:
    currentParams.alpha=rand(1)*0.5+0.25;

    %Initial Value for mu: draw from the normal prior.
    currentParams.mu=options.priors.mu(1)+sqrt(options.priors.mu(2))*randn(1);

    %Initial Value for sigma2: draw from inverse gamma prior is likely a bad
    %ideas, due to the very large possible values. So truncate to below some
    %value:
    while (1)
        t = 1/gamrnd(options.priors.sigma2(1),1/options.priors.sigma2(2));
        if (t < 5)
            currentParams.sigma2 = t;
            break;
        end
    end    

    %Initial Value for phi: draw from the log-normal distrubution, truncated to less
    %than a cutoff, determined by the prior parameters:
    cutt=exp(options.priors.phi(1)+2*sqrt(options.priors.phi(2)));
    while (1)
        t = lognrnd(options.priors.phi(1), sqrt(options.priors.phi(2)));
        if t<cutt
            currentParams.phi=t;
            break;
        end
    end


    %Initial Value for tau2_I: : draw from inverse gamma prior truncated to
    %less than some cut off value:
    while (1)
        t = 1/gamrnd(options.priors.tau2_I(1),1/options.priors.tau2_I(2));
        if t<5
            currentParams.tau2_I=t;
            break;
        end
    end
    
    %Initial Value for tau2_P: : draw from inverse gamma prior truncated to
    %less than some cut off value:
    pars = options.priors.tau2_P;
    vals = zeros(size(pars,1),1);
    for i = 1:size(pars,1)
        while (1)
            t = 1/gamrnd(pars(i,1), 1/pars(i,2));
            if (t < 10)
                vals(i) = t;
                break;
            end
        end
    end
    currentParams.tau2_P = vals;
    
    %Initial value for Beta_1: mode of the normal prior:
    currentParams.Beta_1 = options.priors.Beta_1(:,1) + rand(size(options.priors.Beta_1,1),1).*sqrt(options.priors.Beta_1(:,2));

    %Initial value for Beta_0: mode of normal prior:
    currentParams.Beta_0 = options.priors.Beta_0(:,1) + rand(size(options.priors.Beta_0,1),1).*sqrt(options.priors.Beta_0(:,2));
end

%  Setting the initial values of the temperature matrix.
if options.useModes
    %Method 3: Set all to the mode of the prior for T(0):
    currentField = options.priors.T_0(1)*ones(size(data.instrumental.locations,1),numel(model.timeline));
    
else
    %Method 2: Draw each T(k) from the prior for T(0):
    currentField = options.priors.T_0(1)+sqrt(options.priors.T_0(2))*randn(size(data.instrumental.locations,1), numel(model.timeline));
end