function mu = sampleAutocorrMean(model, params, priors, currentField)
%% Sample autocorrelation mean parameter

postVar = (size(currentField,2)-1)*sum(sum(model.invSpatialCorrMatrix))/params.sigma2;
postVar = 1/(1/priors.mu(2) + (1-params.alpha)^2*postVar);

postMean = sum(currentField(:,2:end)-params.alpha*currentField(:,1:(end-1)),2);
postMean = postVar*((1-params.alpha)*sum(model.invSpatialCorrMatrix*postMean)/params.sigma2 + priors.mu(1)/priors.mu(2));

mu = postMean + sqrt(postVar) * randn;