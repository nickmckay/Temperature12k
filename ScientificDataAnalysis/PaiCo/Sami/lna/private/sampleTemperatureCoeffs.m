function beta = sampleTemperatureCoeffs(data, params, priors, signal)
Cinv = calcInvCovMatrix(numel(signal), params.sigma2L(1), params.phi1L(1), params.phi2L(1));

F = ones(1+numel(data.forcings), numel(signal));
for i = 1:numel(data.forcings)
    F(i+1,:) = data.forcings{i}.data;
end

prioInv = diag(1./priors.beta(:,2));

b = Cinv*F';
A = (F*b + prioInv)^(-0.5);

mn = (signal*b + priors.beta(:,1)'*prioInv)*A*A;

beta = randn(size(mn))*A + mn;
