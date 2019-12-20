function C  = calcCovarianceMatrix(n, sigma2, phi1, phi2)

gamma0 = sigma2 * (1-phi2)/((1+phi2)*(phi1+phi2-1)*(phi2-phi1-1));
rho = zeros(n,1);
rho(1) = 1;
rho(2) = phi1/(1-phi2);
for i = 3:n
    rho(i) = phi1*rho(i-1) + phi2*rho(i-2);
end
C = gamma0*rho(abs(repmat((1:n)',1,n)-repmat(1:n,n,1))+1);