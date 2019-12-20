function field = sampleTempk(data, model, params, k, prevField, nextField)
%% MCMC Sample middle field

m = zeros(size(prevField));
pattern = model.missingPatterns.timePatterns(model.missingPatterns.timeToPattern(k),:);
for j = 1:numel(pattern)
    if (pattern(j)==0)
        % The data does not cover this pattern
        continue; 
    end
    inputPattern = model.missingPatterns.dataPatterns{j}{pattern(j)};
    if (isempty(inputPattern))
        continue;
    end
    if (j == 1)      
        d = data.instrumental.data(inputPattern, data.instrumental.timeIndices == k);
        if (~isempty(d))
            m(inputPattern) = m(inputPattern) + d/params.tau2_I;
        end
    else
        d = data.proxy{j-1}.data(inputPattern.dataIndices, data.proxy{j-1}.timeIndices == k);
        if (~isempty(d))
            m(inputPattern.locationIndices) = m(inputPattern.locationIndices) + ...
                inputPattern.locationMap*(d-params.Beta_0(j-1))*params.Beta_1(j-1)/params.tau2_P(j-1);
        end
    end
end

W = model.invSpatialCorrMatrix/params.sigma2;
m = m + W*(params.alpha*(prevField+nextField) + (1-params.alpha)^2*params.mu);
if (isfield(model,'spatialCovMatrices'))
    sqInvCC = model.sqrtSpatialCovMatrices{model.missingPatterns.timeToPattern(k)};
    invCC = model.spatialCovMatrices{model.missingPatterns.timeToPattern(k)};
else
    [invCC sqInvCC] = calcSingleSpatialCovariance(data, model, params, model.missingPatterns.timeToPattern(k));
end
field = invCC\m + sqInvCC\randn(size(m));

% 
% S = zeros(size(m));
% d = data.instrumental.data(:,data.instrumental.timeIndices == k);
% if (~isempty(d))
%     mask = ~isnan(d);
%     m(mask) = d(mask)/params.tau2_I;
%     S(mask) = 1/params.tau2_I;
% end
% for i = 1:numel(data.proxy)
%     d = data.proxy{i}.data(:,data.proxy{i}.timeIndices == k);
%     mask = ~isnan(d);
%     if (~isempty(d) && any(mask))
%         d = d(mask);
%         ind = data.proxy{i}.locationIndices(mask);
%         for j = 1:numel(ind)
%             m(ind(j)) = m(ind(j)) + ...
%                 (d(j) - params.Beta_0(i))*params.Beta_1(i)/params.tau2_P(i);
%             S(ind(j)) = S(ind(j)) + params.Beta_1(i)^2/params.tau2_P(i);
%         end
%     end
% end
% S = model.invSpatialCorrMatrix*(1+params.alpha^2)/params.sigma2 + diag(S);
% m = m + model.invSpatialCorrMatrix/params.sigma2*(params.alpha*(prevField+nextField) + (1-params.alpha)^2*params.mu);
% field = S\m + chol(S)\randn(size(m));