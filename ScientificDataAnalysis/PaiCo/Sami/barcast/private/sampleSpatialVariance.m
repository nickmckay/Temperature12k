function sigma2 = sampleSpatialVariance(model, params, priors, currentField)
%% Sample autocorrelation mean parameter

postAlpha = (numel(currentField)-size(currentField,1))/2 + priors.sigma2(1);

tDiff = currentField(:,2:end) - params.alpha*currentField(:,1:(end-1)) - (1-params.alpha)*params.mu;
postBeta = priors.sigma2(2) + sum(sum((model.invSpatialCorrMatrix*tDiff).*tDiff))/2;

sigma2 = min(50,1/gamrnd(postAlpha, 1/postBeta));