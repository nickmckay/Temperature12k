function [phi, spatCorr, iSpatCorr] = sampleSpatialCovarianceRange(model, params, priors, mhparams, currentField)
%% Sample spatial covariance range parameter

tDiff = currentField(:,2:end) - params.alpha*currentField(:,1:(end-1)) - (1-params.alpha)*params.mu;

spatCorr = model.spatialCorrMatrix;
logphi = log(params.phi);
value = -1/(2*priors.phi(2))*(logphi - priors.phi(1))^2 ...
    - 1/(2*params.sigma2)*sum(sum(tDiff.*(spatCorr\tDiff)));

propVar = sqrt(mhparams.log_phi(1));
n = (size(currentField,2)-1)/2;

for sample = 1:mhparams.log_phi(2)
    logphiProp = logphi + propVar * randn;
    spatCorrProp = exp(-exp(logphiProp)*model.distances);
    
    valueProp = -1/(2*priors.phi(2))*(logphiProp - priors.phi(1))^2 ...
        - 1/(2*params.sigma2)*sum(sum(tDiff.*(spatCorrProp\tDiff)));
    
    logRatio = valueProp - value + n*log(det(spatCorrProp\spatCorr));

    if (exp(logRatio) > rand)
        % Accept proposition
        logphi = logphiProp;
        spatCorr = spatCorrProp;
        value = valueProp;
    end
end

phi = exp(logphi);
iSpatCorr = inv(spatCorr);