function res = logProbAutocorr(scov, sigma2, phi1, phi2)
% Condensed inverse autocorr matrix
C = [(1+phi1^2+phi2^2);... % Diagonal elements
    (phi1*phi2-phi1); ... % Lag 1 elements
    -phi2; ...  % Lag 2 elements
    -(phi1^2+phi2^2); ... % Fix for boundary elements in 1st and last diagonal positions
    -phi2^2; ... % Fix for boundary elements in 2st and 2nd to last diagonal positions
    -phi1+phi2]; % Fix for boundary elements in ends of off diagonal

res = sum(scov*C)/sigma2;