function field = sampleTemp0(model, params, priors, nextField)
%% MCMC sample values for the field
n = numel(nextField);
C = (params.alpha^2/params.sigma2*model.invSpatialCorrMatrix + eye(n)/priors.T_0(2))^(-1/2);
A = model.invSpatialCorrMatrix*(params.alpha*nextField - params.alpha*(1-params.alpha)*(params.mu*ones(n,1)));
A = A + ones(n,1)/priors.T_0(2)*priors.T_0(1);
A = C*C'*A;
field = A + C*randn(n,1);