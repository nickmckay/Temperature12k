function field = sampleTempField(data, model, params, priors)

nrTime = numel(model.timeline);
nrLoc = size(model.spatialCorrMatrix,1);

% Spatial and temporal correlation

% Matlab removes elements from the top of the spdiags input vector if off
% diagonals are used
AC = spdiags([params.alpha^2;ones(nrTime-2,1)*(1+params.alpha^2);1],0,nrTime,nrTime) + ...    
    spdiags(ones(nrTime,2)*(-params.alpha),[1 -1], nrTime, nrTime); 
S = kron(AC,model.invSpatialCorrMatrix/params.sigma2);
m = (sum(model.invSpatialCorrMatrix,2)/params.sigma2)*...
    [params.alpha*(1-params.alpha) (1-params.alpha)^2*ones(1,nrTime-2) (1-params.alpha)]*params.mu;
m = m(:);

% Add prior for the first time point
inds = 1:nrLoc;
S = S + sparse(inds, inds, ones(numel(inds),1)/priors.T_0(2), size(S,1), size(S,1));
m = m + sparse(inds, ones(numel(inds),1), ones(numel(inds),1)/priors.T_0(2)*priors.T_0(1), size(m,1), 1);

% Add instrumental data
inds = repmat(data.instrumental.locationIndices(:), 1, numel(data.instrumental.timeIndices))+ ...
    repmat(data.instrumental.timeIndices(:)'-1,numel(data.instrumental.locationIndices),1)*nrLoc;
mask = ~isnan(data.instrumental.data);
inds = inds(mask);
S = S + sparse(inds, inds, ones(numel(inds),1)/params.tau2_I, size(S,1), size(S,2)); 
m = m + sparse(inds, ones(numel(inds,1)), data.instrumental.data(mask)/params.tau2_I, size(m,1), 1);

% Add proxy data
for i = 1:numel(data.proxy)
    inds = repmat(data.proxy{i}.locationIndices(:), 1, numel(data.proxy{i}.timeIndices))+ ...
        repmat(data.proxy{i}.timeIndices(:)'-1,numel(data.proxy{i}.locationIndices),1)*nrLoc;
    mask = ~isnan(data.proxy{i}.data);
    inds = inds(mask);

    S = S + sparse(inds, inds, ones(numel(inds),1)*(params.Beta_1(i)^2/params.tau2_P(i)), size(S,1), size(S,2)); 

    % Sparse representation takes care of duplicate indices
    m = m + sparse(inds, ones(numel(inds,1)), ...
        (data.proxy{i}.data(mask)-params.Beta_0(i))*params.Beta_1(i)/params.tau2_P(i), size(m,1), 1);
end

% Sample field
field = S\m + chol(S)\randn(size(m));
field = reshape(field,nrLoc, nrTime);