function [beta0 beta1] = sampleProxyCoeff(data, params, priors, currentField, proxyId)
%% Sample proxy measurement multiplier

field = currentField(data.proxy{proxyId}.locationIndices, data.proxy{proxyId}.timeIndices);
mask = ~isnan(data.proxy{proxyId}.data);
field = field(mask);
field = field(:);
data = data.proxy{proxyId}.data(mask);
data = data(:);

% Prior
S = [1/priors.Beta_0(proxyId, 2) 0; 0 1/priors.Beta_1(proxyId,2)];
m = S*[priors.Beta_0(proxyId, 1); priors.Beta_1(proxyId,1)];

% Data
X = [ones(size(field)) field];
S = S + X'*X/params.tau2_P(proxyId);
m = m + X'*data/params.tau2_P(proxyId);

coeff = S\m + chol(S)\randn(size(m));
beta0 = coeff(1);
beta1 = coeff(2);


