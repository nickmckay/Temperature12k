function tau2 = sampleInstrumentalErrorVar(data, priors, currentField)
%% Sample instrumental measurement error variance

res = data.instrumental.data - currentField(data.instrumental.locationIndices, data.instrumental.timeIndices);
res = res(~isnan(res));
res = res(:);
postBeta = res'*res/2 + priors.tau2_I(2);
% As currentField is completely filled with no NaN's, the size of res
% equals the number of values in instrumental data
postAlpha = numel(res)/2 + priors.tau2_I(1);

tau2 = min(50,1/gamrnd(postAlpha, 1/postBeta));