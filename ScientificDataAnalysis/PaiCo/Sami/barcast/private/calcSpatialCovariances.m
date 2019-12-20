function [covMtxs, sqrtCovMtxs] = calcSpatialCovariances(data, model, params)
%% Calculates the temporal covariance matries
% Done for each pattern of missing data

covMtxs = cell(size(model.missingPatterns.timePatterns,1),1);
sqrtCovMtxs = covMtxs;
n = size(data.instrumental.locations,1);
for i = 1:size(model.missingPatterns.timePatterns,1)
    C = zeros(n,1);
    for j = 1:size(model.missingPatterns.timePatterns,2)
        if (model.missingPatterns.timePatterns(i,j)==0)
            % The data does not cover this pattern
            continue; 
        end
        pattern = model.missingPatterns.dataPatterns{j}{model.missingPatterns.timePatterns(i,j)};
        if (j == 1)
            % For instrumental data, pattern contains the indices of the
            % locations that have data for timePattern            
            C(pattern) = C(pattern) + 1/params.tau2_I;
        else            
            C(pattern.locationIndices) = C(pattern.locationIndices) + ...
                sum(pattern.locationMap,2)*params.Beta_1(j-1)^2/params.tau2_P(j-1);
        end        
    end
    
    C = model.invSpatialCorrMatrix*(1+params.alpha^2)/params.sigma2 + diag(C);
%     fprintf('.'); drawnow;
%     disp([num2str(i) '/' num2str(size(model.missingPatterns.timePatterns,1))]);
%     memory
    
    sqrtCovMtxs{i} = chol(C);
    covMtxs{i} = C;
end
% fprintf('\n');