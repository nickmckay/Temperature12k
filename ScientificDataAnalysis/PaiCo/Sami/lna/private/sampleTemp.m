function signal = sampleTemp(data, model, params)

% % The code in Li et al. 2009 was missing distribution involving forcings
% % from the mean of the temperature distribution
%
% Cinv = calcInvCovMatrix(numel(model.timeline), params.sigma2L(1), params.phi1L(1), params.phi2L(1));
% Ainv = Cinv(model.missingMask, model.missingMask);
% b = zeros(1,size(Ainv,1));
% 
% for i = 1:numel(data.proxy)
%     Cinv = calcInvCovMatrix(size(data.proxy{i}.data,2), params.sigma2L(1+i), params.phi1L(1+i), params.phi2L(1+i));
%     G = data.proxy{i}.transform*Cinv;
%     H = G*data.proxy{i}.transform';
%     
%     Ainv = Ainv + sum(params.betaL{i}.^2)*H(model.missingMask, model.missingMask);
%     b = b + sum((data.proxy{i}.data-repmat(params.muL{i},1,size(data.proxy{i}.data,2))).* ...
%         repmat(params.betaL{i},1,size(data.proxy{i}.data,2))*G(model.missingMask,:)',1);
%     b = b - sum(params.betaL{i}.^2)*model.instrumental*H(~model.missingMask,model.missingMask);
% end
% 
% signal = zeros(1,numel(model.timeline));
% A = Ainv^(-0.5);
% signal(model.missingMask) = b*A*A+ randn(1,nnz(model.missingMask))*A;
% signal(data.instrumental.timeIndices) = model.instrumental;

% Fixed version

% Instrumental data
Qinv = calcInvCovMatrix(numel(data.target.times), params.sigma2L(1), params.phi1L(1), params.phi2L(1));
muhat = sum(Qinv,2)*params.beta(1);        
% for j = 1:numel(data.forcings)
%     muhat = muhat + params.beta(j+1)*data.forcings{j}.data;
% end

for i = 1:numel(data.proxy)
    Cj = calcInvCovMatrix(size(data.proxy{i}.data,2), params.sigma2L(1+i), params.phi1L(1+i), params.phi2L(1+i));

    Cj = data.proxy{i}.transform*Cj;


    Qinv = Qinv + sum(params.betaL{i}.^2)*Cj*data.proxy{i}.transform';
    muhat = muhat + Cj*(data.proxy{i}.data - repmat(params.muL{i},1,size(data.proxy{i}.data,2)))'*params.betaL{i};
end

S = Qinv(model.missingMask, model.missingMask);
m = muhat(model.missingMask) - Qinv(model.missingMask,~model.missingMask)*model.instrumental';

signal = zeros(1, numel(data.target.times));
signal(~model.missingMask) = model.instrumental;
signal(model.missingMask) = S\m + chol(S)\randn(size(m));