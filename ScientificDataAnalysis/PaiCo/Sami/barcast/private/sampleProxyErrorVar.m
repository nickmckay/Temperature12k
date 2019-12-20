function tau2 = sampleProxyErrorVar(data, params, priors, currentField, proxyId)
%% Sample instrumental measurement error variance

res = params.Beta_1(proxyId)*currentField(data.proxy{proxyId}.locationIndices, data.proxy{proxyId}.timeIndices);
res = data.proxy{proxyId}.data - (res + params.Beta_0(proxyId));
res = res(~isnan(res));
res = res(:);
postBeta = res'*res/2 + priors.tau2_P(proxyId,2);
% As currentField is completely filled with numbers, the size of res
% equals the number of non-NaN's in instrumental data
postAlpha = numel(res)/2 + priors.tau2_P(proxyId,1);

tau2 = min(5,1/gamrnd(postAlpha, 1/postBeta));