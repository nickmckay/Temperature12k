function [B] = nanzscore(A)

B = (A-nanmean(A))./nanstd(A);