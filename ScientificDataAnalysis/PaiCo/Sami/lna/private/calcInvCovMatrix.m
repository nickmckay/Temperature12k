function Cinv = calcInvCovMatrix(n, sigma2, phi1, phi2)

% Cinv = diag(ones(n,1)*(1+phi1^2+phi2^2));
% Cinv = Cinv + (phi1*phi2-phi1)*(diag(ones(n-1,1),1)+diag(ones(n-1,1),-1));
% Cinv = Cinv - phi2*(diag(ones(n-2,1),2)+diag(ones(n-2,1),-2));
% Cinv(1:2,1:2) = [1, -phi1; -phi1, 1+phi1^2];
% Cinv((end-1):end, (end-1):end) = [1+phi1^2, -phi1; -phi1, 1];
% CO = Cinv / sigma2;

% Using sparse representation is faster

Cinv = full(spdiags([...
    ones(n,1)*(1+phi1^2+phi2^2), ... % Diagonal
    ones(n,2)*(phi1*phi2-phi1), ... % 1 off diaginal
    ones(n,2)*-phi2]/sigma2, ... % 2 off diagonal
    [0 1 -1 2 -2], n, n));

Cinv(1:2,1:2) = [1, -phi1; -phi1, 1+phi1^2]/sigma2;
Cinv((end-1):end, (end-1):end) = [1+phi1^2, -phi1; -phi1, 1]/sigma2;