function beta1 = sampleProxyMultiplier(data, params, priors, currentField, proxyId)
%% Sample proxy measurement multiplier

field = currentField(data.proxy{proxyId}.locationIndices, data.proxy{proxyId}.timeIndices);
mask = ~isnan(data.proxy{proxyId}.data);
field = field(mask);
field = field(:);
data = data.proxy{proxyId}.data(mask);
data = data(:);

postVar = 1/(1/priors.Beta_1(proxyId,2) + field'*field/params.tau2_P(proxyId));
postMean = postVar * (priors.Beta_1(proxyId,1)/priors.Beta_1(proxyId,2) ...
    + (data-params.Beta_0(proxyId))'*field/params.tau2_P(proxyId));

beta1 = postMean + sqrt(postVar)*randn;