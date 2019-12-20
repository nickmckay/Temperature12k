function [invCC sqInvCC] = calcSingleSpatialCovariance(data, model, params, timePattern)

pattern = model.missingPatterns.timePatterns(timePattern,:);
C = zeros(size(data.instrumental.locations,1),1);
for j = 1:numel(pattern)
    if (pattern(j)==0)
        % The data does not cover this pattern
        continue; 
    end
    inputPattern = model.missingPatterns.dataPatterns{j}{pattern(j)};
    if (j == 1)
        % For instrumental data, pattern contains the indices of the
        % locations that have data for timePattern
        C(inputPattern) = C(inputPattern) + 1/params.tau2_I;
    else
        C(inputPattern.locationIndices) = C(inputPattern.locationIndices) + ...
            sum(inputPattern.locationMap,2)*params.Beta_1(j-1)^2/params.tau2_P(j-1);
    end        
end

invCC = model.invSpatialCorrMatrix*(1+params.alpha^2)/params.sigma2 + diag(C);
sqInvCC = chol(invCC);