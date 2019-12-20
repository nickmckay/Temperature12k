function sigma2 = sampleSigma2(data, params, priors, signal)

sigma2 = params.sigma2L;


for i = 1:(numel(data.proxy)+1)

    if (i == 1)
        [m,n] = size(signal);
        Cinv = calcInvCovMatrix(n, 1, params.phi1L(i), params.phi2L(i));    
        R = signal - params.beta(1);
        for j = 1:numel(data.forcings)
            R = R - params.beta(j+1)*data.forcings{j}.data;
        end
        r = R*Cinv*R';
    else

        
        [m,n] = size(data.proxy{i-1}.data);
        Cinv = calcInvCovMatrix(n, 1, params.phi1L(i), params.phi2L(i));    
        R = data.proxy{i-1}.data - repmat(params.muL{i-1},1,n);
        R = R - params.betaL{i-1}*(signal*data.proxy{i-1}.transform);
        
        r = sum(sum((R*Cinv).*R));
%         r = trace(Cinv*(R'*R));
    end
    
    
    alpha = priors.sigma2L(i,1) + n*m/2;
    rhat = r/2 + priors.sigma2L(i,2)^(-1);
    
    sigma2(i) = 1./gamrnd(alpha, 1./rhat);
end