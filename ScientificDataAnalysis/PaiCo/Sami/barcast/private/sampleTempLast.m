function field = sampleTempLast(data, model, params, prevField)
%% MCMC Sample last field

field = zeros(size(prevField));
pattern = model.missingPatterns.timePatterns(model.missingPatterns.timeToPattern(end),:);
C = zeros(numel(prevField),1);
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
        % Update spatial covariance matrix temp variable
        C(inputPattern) = C(inputPattern) + 1/params.tau2_I;

        % Update post mean helper variable
        d = data.instrumental.data(inputPattern, data.instrumental.timeIndices == model.timeline(end));
        if (~isempty(d))
            field(inputPattern) = field(inputPattern) + d/params.tau2_I;
        end
    else
        % Update spatial covariance matrix temp variable
        C(inputPattern.locationIndices) = C(inputPattern.locationIndices) + ...
            sum(inputPattern.locationMap,2)*params.Beta_1(j-1)^2/params.tau2_P(j-1);

        % Update post mean helper variable
        d = data.proxy{j-1}.data(inputPattern.dataIndices, data.proxy{j-1}.timeIndices == model.timeline(end));
        if (~isempty(d))
            field(inputPattern.locationIndices) = field(inputPattern.locationIndices) + ...
                inputPattern.locationMap*params.Beta_1(j-1)/params.tau2_P(j-1)*(d-params.Beta_0(j-1));
        end
    end
 end

S = model.invSpatialCorrMatrix/params.sigma2 + diag(C);
m = model.invSpatialCorrMatrix/params.sigma2*(params.alpha*prevField + (1-params.alpha)*params.mu);
field = S\(field + m) + chol(S)\randn(size(field));