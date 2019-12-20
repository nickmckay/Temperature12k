function [phi1, phi2] = sampleAutocorr(data, params, mhparams, signal)

phi1 = params.phi1L;
phi2 = params.phi2L;

propStd = sqrt(mhparams.phi(1));
lapcdf = @(x,mu)(.5*(1+sign(x-mu)*(1-exp(-abs(x-mu)/propStd))));
lapinv = @(p,mu)(mu - propStd*sign(p-.5)*log(1-2*abs(p-.5)));
lappdf = @(x,mu)(exp(-abs(x-mu)/propStd)/(2*propStd));

for i = 1:(numel(data.proxy)+1)

    p1 = phi1(i);
    p2 = phi2(i);
        
    if (i == 1)
        R = signal - params.beta(1);        
        for j = 1:numel(data.forcings)
            R = R - params.beta(j+1)*data.forcings{j}.data;
        end
        [m,n] = size(signal);
    else
        [m,n] = size(data.proxy{i-1}.data);        
        R = data.proxy{i-1}.data - repmat(params.muL{i-1},1,n);
        R = R - params.betaL{i-1}*(signal*data.proxy{i-1}.transform);        
    end
    
%     Cinv = calcInvCovMatrix(n, params.sigma2L(i), p1, p2);
    logdet = calcLogDetInvCovMatrix(n, params.sigma2L(i), p1, p2);
   
    
%     store = zeros(3,mhparams.phi(2));    
    
    sampleCov = zeros(m,6);
    % Lag 0
    sampleCov(:,1) = sum(R.*R,2);
    % Lag 1
    sampleCov(:,2) = 2*sum(R(:,2:end).*R(:,1:(end-1)),2);
    % Lag 2
    sampleCov(:,3) = 2*sum(R(:,3:end).*R(:,1:(end-2)),2);    
    % Helpers to handle boundaries
    sampleCov(:,4) = R(:,1).^2 + R(:,end).^2;
    sampleCov(:,5) = R(:,2).^2 + R(:,end-1).^2;
    sampleCov(:,6) = R(:,1).*R(:,2) + R(:,(end-1)).*R(:,end);

    res = logProbAutocorr(sampleCov, params.sigma2L(i), p1, p2);
    
%     % Plot log pdf in a grid
%     [XX, YY] = meshgrid(linspace(-1,1,1e2), linspace(-1,1,1e2));
%     ZZ = zeros(size(XX));
%     for x = 1:size(XX,1), 
%         for y = 1:size(XX,2), 
%             ld = calcLogDetInvCovMatrix(n, params.sigma2L(i), XX(x,y), YY(x,y)); 
%             nr =  logProbAutocorr(sampleCov, params.sigma2L(i), XX(x,y), YY(x,y)); 
%             ZZ(x,y) = m*ld - nr; 
%         end
%     end
%     ZZ(imag(ZZ) > 0) = nan;
%     ZZ(abs(ZZ) == inf) = nan;
%     ZZ = real(ZZ);
%     figure(3+i); clf; 
%     surf(XX,YY,exp(ZZ-max(max(ZZ))));
%     pause;

    for sample = 1:mhparams.phi(2)
%         newp2 = p2 + randn*propStd;
%         while (abs(newp2) > 1)
%             newp2 = p2 + randn*propStd;
%         end
%         
%         if (abs(newp2) < 0.99)
%             newp1 = p1 + randn*propStd;
%             while (abs(newp1) > 1-abs(newp2))
%                 newp1 = p1 + randn*propStd;
%             end
%         else
%             newp1 = (2*rand - 1)*(1-abs(newp2));
%         end

        % Laplace distribution is easier to integrate
        newp2 = lapinv(rand*(lapcdf(1,p2) - lapcdf(-1,p2)) + lapcdf(-1,p2),p2);        
        newp1 = lapinv(rand*(lapcdf(1-abs(newp2),p1) - lapcdf(-1+abs(newp2), p1)) + lapcdf(-1+abs(newp2), p1), p1);
        
%         newp2 = rand*(limp2(2)-limp2(1)) + limp2(1);
%         
%         if (abs(newp2) < 0.99)
%             newp1 = p1 + (rand-.5)*propStd;
%             while (abs(newp1) > 1-abs(newp2))
%                 newp1 = p1 + (rand-.5)*propStd;
%             end
%         else
%             newp1 = (2*rand - 1)*(1-abs(newp2));
%         end
%         
%         if (1-abs(newp2) > 1e-9)
%             prob = [.5*erfc(-(-1+abs(newp2)-p1)/propStd/sqrt(2)) .5*erfc(-(1-abs(newp2) - p1)/propStd/sqrt(2))]; 
% %             prob = [normcdf(-1+abs(newp2), p1, propStd) normcdf(1-abs(newp2), p1, propStd)];
%             if (prob(2)-prob(1) > 1e-9)
%                 newp1 = -sqrt(2)*erfcinv(2*(rand*(prob(2)-prob(1)) + prob(1)))*propStd + p2;                
% %                 newp1 = norminv(rand*(prob(2)-prob(1)) + prob(1), p1, propStd);
%             else
%                 newp1 = (2*rand-1)*(1-abs(newp2));
%             end
%         else
%             newp1 = (2*rand-1)*(1-abs(newp2));
%         end
%         

        newlogdet = calcLogDetInvCovMatrix(n, params.sigma2L(i), newp1, newp2);    
        newres = logProbAutocorr(sampleCov, params.sigma2L(i), newp1, newp2);        
                
        % Li et al. 2009 log accept probability. The probabilities for
        % proposal distribution are incorrect.
%         Cnew = calcInvCovMatrix(n, params.sigma2L(i), newp1, newp2);
%         logr1 = -m*log(det(Cnew)) + R*Cnew*R' + ... 
%                 log((1-p2)*(normcdf(1,newp2, propStd) - normcdf(-1, newp2, propStd)) + ...
%                 abs(normcdf(1-p2, newp1, propStd)-normcdf(p2-1, newp1, propStd)));
%             
%         logr2 = -m*log(det(Cinv)) + R*Cinv*R' + ... 
%                 log((1-newp2)*(normcdf(1, p2, propStd) - normcdf(-1, p2, propStd)) + ...
%                 abs(normcdf(1-newp2, p1, propStd) - normcdf(newp2 - 1, p1, propStd)));
%         Q = log((1-p2)*(normcdf(1,newp2, propStd) - normcdf(-1, newp2, propStd)) + ...
%                 abs(normcdf(1-p2, newp1, propStd)-normcdf(p2-1, newp1, propStd))) + ...
%                 log((1-newp2)*(normcdf(1, p2, propStd) - normcdf(-1, p2, propStd)) + ...
%                 abs(normcdf(1-newp2, p1, propStd) - normcdf(newp2 - 1, p1, propStd)));

        % Fixed version                
        
        % PDF ratio
        logr = 0.5*m*(newlogdet - logdet);
        logr = logr - 0.5*(newres-res);
        
%         % Log proposal density from newp2 to p2
%         Qnp2p2 = log(normpdf(p2, newp2, propStd)) - log(normcdf(1, newp2, propStd) - normcdf(-1, newp2, propStd));
%         Qnp2p2 = log(normpdf(p2, newp2, propStd)) - log(.5*erfc(-(1-newp2)/propStd/sqrt(2)) - .5*erfc(-(-1-newp2)/propStd/sqrt(2)));
%         Qnp2p2 = log(min(1, newp2 + propStd/2) - max(0, newp2 - propStd/2));          
        Qnp2p2 = log(lappdf(p2, newp2)) - log(lapcdf(1, newp2) - lapcdf(-1, newp2));
        
        % Log proposal density from p2 to newp2
%         Qp2np2 = log(normpdf(newp2, p2, propStd)) - log(normcdf(1, p2, propStd) - normcdf(-1, p2, propStd));
%         Qp2np2 = log(normpdf(newp2, p2, propStd)) - log(.5*erfc(-(1-p2)/propStd/sqrt(2)) - .5*erfc(-(-1-p2)/propStd/sqrt(2)));
%         Qp2np2 = log(min(1, p2 + propStd/2) - max(0, p2 - propStd/2));
        Qp2np2 = log(lappdf(newp2, p2)) - log(lapcdf(1, p2) - lapcdf(-1, p2));

        % Log proposal density from newp1 to p1 given newp2
%         Qnp1p1 = log(normpdf(p1, newp1, propStd)) - log(normcdf(1-abs(newp2), newp1, propStd) - normcdf(-1+abs(newp2), newp1, propStd));
%         Qnp1p1 = log(normpdf(p1, newp1, propStd)) - log(.5*erfc(-(1-abs(newp2)-newp1)/propStd/sqrt(2)) - .5*erfc(-(-1+abs(newp2)-newp1)/propStd/sqrt(2)));
%         Qnp1p1 = log(min(1-abs(newp2), newp1 + propStd/2) - max(-1+abs(newp2), newp1 - propStd/2));
        Qnp1p1 = log(lappdf(p1, newp1)) - log(lapcdf(1-abs(newp2), newp1) - lapcdf(-1+abs(newp2), newp1));

        % Log proposal density from p1 to newp1 given newp2
%         Qp1np1 = log(normpdf(newp1, p1, propStd)) - log(normcdf(1-abs(newp2), p1, propStd) - normcdf(-1+abs(newp2), p1, propStd)); 
%         Qp1np1 = log(normpdf(newp1, p1, propStd)) - log(.5*erfc(-(1-abs(newp2)-p1)/propStd/sqrt(2)) - .5*erfc(-(-1+abs(newp2)-p1)/propStd/sqrt(2)));
%         Qp1np1 = log(min(1-abs(newp2), p1 + propStd/2) - max(-1+abs(newp2), p1 - propStd/2));
        Qp1np1 = log(lappdf(newp1, p1)) - log(lapcdf(1-abs(newp2), p1) - lapcdf(-1+abs(newp2), p1));
        
        logr = logr + Qnp2p2 - Qp2np2 + Qnp1p1 - Qp1np1;

        
        if (log(rand) < logr)
            % Accept new phi's
            p1 = newp1;
            p2 = newp2;
%             Cinv = Cnew;
            logdet = newlogdet;
            res = newres;
        end        
    end
    
%     % Comparison to Yule-Walker ML cofficients
%     if (i == 1)
%         p = polydata(ar(R,2));
%         disp([num2str(i) ': p1=' num2str(p1) ' arp1=' num2str(-p(2)) ' p2=' num2str(p2) ' arp2=' num2str(-p(3))]);
%     end
    
    phi1(i) = p1;
    phi2(i) = p2;
end