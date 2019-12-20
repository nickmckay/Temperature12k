function mbest = eigenselect(eigen)
%    EIGENSELECT(series) finds optimal choice of # of eigenvalues
%    to retain based on estimate of noise background trend to
%    log eigenvalue spectrum
%
%    restrict to space spanned by first 0.5*M eigenvalues (to avoid
%    non-linear regime at high eigenvalue number)
%
%    From Mann et al. 2007


eigen = sort(eigen/sum(eigen),'descend');

M = round(numel(eigen)/2);
A = [ones(M,1) (1:M)'];
series = log(eigen(1:M));
beta = A\series;
trend = A*beta;
%    Keep the set of first K eigenvalues where Kth eigenvalue
%    is the first to drop below background trend within indicated
%    tolerance
%
%    use tolerance of 1 standard deviation of trend residual
%   
mbest = find(series < (trend + std(trend-series)),1,'first');