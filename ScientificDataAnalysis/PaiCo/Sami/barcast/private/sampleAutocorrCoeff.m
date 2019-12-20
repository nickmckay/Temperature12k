function alpha = sampleAutocorrCoeff(model, params, priors, currentField)
%% Sample autocorrelation coefficient

currentField = currentField - params.mu;
field = currentField(:,1:(end-1))'*model.invSpatialCorrMatrix/params.sigma2;
postVar = 1/sum(sum(field.*currentField(:,1:(end-1))'));
postMean = postVar * sum(sum(field.*currentField(:,2:end)'));
postStd = sqrt(postVar);

prob = normcdf(priors.alpha, postMean, postStd);
if (prob(2)-prob(1) < 1e-9)
    alpha = rand*(priors.alpha(2)-priors.alpha(1)) + priors.alpha(1);
else
    alpha = norminv((prob(2)-prob(1))*rand + prob(1), postMean, postStd);
end
% while 1
%     alpha = postMean + postStd*randn;
%     if (alpha > priors.alpha(1) && alpha < priors.alpha(2))
%         break;
%     end
% end