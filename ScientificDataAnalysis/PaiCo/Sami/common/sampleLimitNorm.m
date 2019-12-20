function x = sampleLimitNorm(m, s, limits)

cdfL = normcdf(limits(1),m,s);
cdfR = normcdf(limits(2),m,s);

if ((cdfL == 0 && cdfR == 0) || (cdfL == 1 && cdfR == 1))
    x = m;
else
    x = norminv(rand*(cdfR-cdfL)+cdfL,m,s);
end

if (x < limits(1))
    x = limits(1);
elseif (x > limits(2))
    x = limits(2);
end
