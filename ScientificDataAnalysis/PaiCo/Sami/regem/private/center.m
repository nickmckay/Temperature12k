function [X, M] = center(X)

M = nanmean(X,1);
X = X - repmat(M,size(X,1),1);