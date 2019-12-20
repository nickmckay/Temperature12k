function [muL, betaL] = sampleProxyCoeffs(data, params, priors, signal)

muL = params.muL;
betaL = params.betaL;

for i = 1:numel(data.proxy)
    Cinv = calcInvCovMatrix(size(data.proxy{i}.data,2), params.sigma2L(1+i), params.phi1L(1+i), params.phi2L(1+i));
    
    MT = [ones(1,size(data.proxy{i}.data,2));signal*data.proxy{i}.transform];
    prioInv = diag(1./[priors.muL(i,2) priors.betaL(i,2)]);
    
    b = Cinv*MT';
    A = MT*b + prioInv;
   
    mn = (data.proxy{i}.data*b)' + repmat(prioInv*[priors.muL(i,1);priors.betaL(i,1)],1,size(data.proxy{i}.data,1));
    
    r = A\mn + chol(A)\randn(size(mn));
    muL{i} = r(1,:)';
    betaL{i} = r(2,:)';
end