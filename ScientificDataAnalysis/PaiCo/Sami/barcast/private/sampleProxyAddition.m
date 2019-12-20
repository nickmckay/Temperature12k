function beta0 = sampleProxyAddition(data, params, priors, currentField, proxyId)
%% Sample proxy measurement addition parameter

res = data.proxy{proxyId}.data ...
    - params.Beta_1(proxyId)*currentField(data.proxy{proxyId}.locationIndices, data.proxy{proxyId}.timeIndices);
res = res(~isnan(res));

postVar = 1/(1/priors.Beta_0(proxyId,2) + numel(res)/params.tau2_P(proxyId));
postMean = postVar*(priors.Beta_0(proxyId,1)/priors.Beta_0(proxyId,2) ...
    + sum(res)/params.tau2_P(proxyId));

beta0 = postMean + sqrt(postVar) * randn;