function d = calcLogDetInvCovMatrix(n, sigma2, phi1, phi2)
d = -n*log(sigma2);
d = d + log((1-phi2^2)^2 - phi1^2*(1+phi2)^2);