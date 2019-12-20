function [bestinitf bestinitlogl] = optimize(A, counts, options)
% Newton-Raphson iteration to find maximum likelihood estimate for f

if (options.heuristicStart)
    initf = full(((counts(:,1)-counts(:,2))./sum(counts,2))'*A)';
    initf = {initf/std(initf)};
else       
    initf = cell(options.randomStarts,1);
    for i = 1:options.randomStarts
        initf{i} = rand(size(A,2),1);
    end
end

if (numel(initf) > 1)
    tracker('PAICO>randomStarts');
end

bestinitlogl = -inf;
bestinitf = initf{1};

for bi = 1:numel(initf)

    f = initf{bi};
    oldlogl = NaN;
    logl = -inf;
    iters = 0;
    errstd = 1;
    z = A*f;

    bestlogl = -inf;
    bestf = f;
    cdf = normcdf(z,0,errstd);
    cdf(cdf==0) = eps;
    cdf(cdf==1) = 1-eps;
    if (isempty(options.damping))
        v = 0;
    else
        v = options.damping;
    end
    
    while ~(abs(logl-oldlogl)/abs(oldlogl) < options.errorTolerance) && (iters < options.maxIters) && ~isnan(logl)%  
        oldlogl = logl;

        % Newton-Raphson method to find maximum likelihood
        pdf = normpdf(z,0,errstd);
        
        d = counts(:,1).*(pdf./cdf.*z + (pdf./cdf).^2);
        d = d + counts(:,2).*(pdf./(1-cdf).*(-z) + (pdf./(1-cdf)).^2);
        H = -A'*sparse(1:numel(z),1:numel(z),d,numel(z),numel(z))*A - options.regcov;
        
        F = ((counts(:,1).*(pdf./cdf) - counts(:,2).*(pdf./(1-cdf)))'*A)';
        
        % Levenberg-Marquardt style damping if X'*X is close to singular.
        f = f - (H+diag(diag(H))*v)\(F - options.regcov*f);

        z = A*f;
        cdf = normcdf(z,0,errstd);
        % Fix numerical instability issues
        cdf(cdf==0) = eps;
        cdf(cdf==1) = 1-eps;
        logl = counts(:,1)'*log(cdf)+counts(:,2)'*log(1-cdf);

        if (logl > bestlogl)
            bestlogl = logl;
            bestf = f;            
        else
            % Adaptive damping increases the damping scalar when
            % Newton-Raphson would decrease the log likelihood. Slows down
            % the convergence, but makes it more stable.             
            v = max(v*10,1e-2);
%             fprintf('PAICO>Optimize: Increased damping to %d\n',v);
            f = bestf;
            logl = bestlogl;
            oldlogl = logl*1.1;
            z = A*f;
            cdf = normcdf(z,0,errstd);
            % Fix numerical instability issues
            cdf(cdf==0) = eps;
            cdf(cdf==1) = 1-eps;            
        end
        
%         figure(1); clf; 
%         plot(f);

        iters = iters + 1;

%         fprintf('PAICO>Optimize: logl=%5.15f diff=%5.15f iter=%d stdb=%5.9f\n', logl, (logl-oldlogl)/abs(oldlogl),iters, std(f));
    end

    if (bestlogl > bestinitlogl)
        bestinitlogl = bestlogl;
        bestinitf = bestf;
    end
    
    if (numel(initf) > 1)
        tracker('PAICO>randomStarts', bi, numel(initf));
    end    
end