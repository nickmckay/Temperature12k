function result = lna(data, options, tracker)
%LNA Multiproxy climate reconstruction with Bayesian Hierarchical models
%   result = lna(data, options, tracker) calculates samples of the climate 
%   signals
%
%   If options is omitted, defaults will be used. If tracker is omitted, no
%   output is generated.
%   
%   If data is set to 'info', result returns default options and
%   information about their meaning.
%
%   For information about data and result structures, type:
%   help rtstruct 
%
%   Method is based on article:
%   Li B, Nychka DW, Ammann CM (2010) The value of multiproxy 
%   reconstruction of past climate. Journal of the American Statistical
%   Association  105(491):883–895, DOI 10.1198/jasa.2010.ap09379.
%
%   Ported to Matlab by Sami Hanhijärvi (2011)


%   COPYRIGHT 
%       2012, Sami Hanhijärvi <sami.hanhijarvi@iki.fi>
%
%   LICENSE
%       This program is free software: you can redistribute it and/or modify
%       it under the terms of the GNU General Public License as published by
%       the Free Software Foundation, either version 3 of the License, or
%       (at your option) any later version.
%
%       This program is distributed in the hope that it will be useful,
%       but WITHOUT ANY WARRANTY; without even the implied warranty of
%       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%       GNU General Public License for more details.
%
%       You should have received a copy of the GNU General Public License
%       along with this program.  If not, see <http://www.gnu.org/licenses/>.

if nargin == 0
    error('RECON:LNA','At least one input parameter has to be defined.');    
end

if (ischar(data) && strcmp(data,'info'))
    % User requested
    [options info] = defaults([]);
    result.options = options;
    result.info = info;
    return;
end

% Check data
rtstruct(data);

% Check options and fill in values that are not defined
if (~exist('options','var'))
    options = [];
end
options = defaults(data, options);  

if(~exist('tracker','var'))
    % Empty function
    tracker = @(a,b,c)([]);
end

%% Initialize
if (~isfield(data, 'forcings'))
    data.forcings = [];
end

model.timeline = data.target.times;

for i = 1:numel(data.proxy)
    % Check that the proxy contains no NaN's as we can not handle them
    mask = ~any(isnan(data.proxy{i}.data),1);
    if (any(~mask))
        data.proxy{i}.data = data.proxy{i}.data(:,mask);
        data.proxy{i}.times = data.proxy{i}.times(mask);
        if (isfield(data.proxy{i},'transform'))
            data.proxy{i}.transform = data.proxy{i}.transform(:,mask);
        end
    end    
    
    % Find which times in proxies correspond to timeline points
    [mask, loc] = ismember(model.timeline(:), data.proxy{i}.times(:));
%     data.proxy{i}.timeIndices = loc;
         
    % Check also that proxies have transformations defined
    % NOTICE: We use proxies as row vectors, whileas Li et al. use column
    % vectors. Therefore, all vectors are transformed in this context.
    if (~isfield(data.proxy{i},'transform'))
        inTargetRes = false;
        if (isfield(data.proxy{i},'intargetres'))
            inTargetRes = data.proxy{pi}.intargetres;
        end

        if (~isfield(data.proxy{i},'lower'))
            % Proxy record doesn't contain time frames for samples. Assume it
            % is annual.
            inTargetRes = true;
        end
        
        if (inTargetRes)

            T = zeros(numel(model.timeline), size(data.proxy{i}.data,2));
            T(sub2ind(size(T),find(mask),loc(loc>0))) = 1;
            data.proxy{i}.transform = T;
        else
            T = zeros(numel(model.timeline), size(data.proxy{i}.data,2));
            
            % Fractions of years are rounded to closest integer year. Higher
            % resolutions than the target can not be reasonably handled with 
            % this implementation.
            lower = round(data.proxy{i}.lower);
            upper = round(data.proxy{i}.upper);

            % Sort data so that times are ascending
%             [lower,ind] = sort(lower,'ascend');
%             upper = upper(ind);

            for j = 1:numel(lower)
                % Lower bound in time is inclusive, upper bound is exclusive.
                mask = model.timeline >= lower(j) & model.timeline < upper(j);
                T(mask,j) = 1/nnz(mask);
            end
            data.proxy{i}.transform = T;
        end
        
    end    
    
    % If the transformations are sparse (<15% nonzero), produce sparse 
    % matrices out of them. 15% was found to be the pivotal value with 
    % simple empirical trial resempling calculations in sampleTemp.
    if (nnz(data.proxy{i}.transform) < .15*numel(data.proxy{i}.transform))
        data.proxy{i}.transform = sparse(data.proxy{i}.transform);
    end
end
% data.instrumental.timeIndices = ismember(model.timeline(:), data.instrumental.times(:));
model.instrumental = nanmean(data.instrumental.data,1);
mask = ismember(data.instrumental.times, model.timeline);
mask = mask & ~isnan(model.instrumental);
model.instrumental = model.instrumental(mask);
model.missingMask = ~ismember(model.timeline, data.instrumental.times(mask));

% Draw initial values for the Bayesian model
[currentParams, currentSignal] = initialValues(data, model, options);

tracker('LNA');
totalIterations = options.preSamplerIterations + options.samplerIterations;

signals = zeros(options.samplerIterations, numel(currentSignal));
presignals = zeros(options.preSamplerIterations, numel(currentSignal));
params = currentParams;
params(totalIterations).beta = [];

%% Carry out actual sampling
for sample = 1:totalIterations
    
    % Sample temperature signal
    currentSignal = sampleTemp(data, model, currentParams);       

    % Sample proxy transformation coefficients
    [currentParams.muL currentParams.betaL] = sampleProxyCoeffs(data, currentParams, options.priors, currentSignal);
    
    % Sample temperature and forcings coefficients
    currentParams.beta = sampleTemperatureCoeffs(data, currentParams, options.priors, currentSignal);            
    
    % Sample autocorrelation coefficients
    [currentParams.phi1L, currentParams.phi2L] = sampleAutocorr(data, currentParams, options.MHparams, currentSignal);
    
    % Sample residual variances
    currentParams.sigma2L = sampleSigma2(data, currentParams, options.priors, currentSignal);
    
    params(sample) = currentParams;    
    
    if (sample > options.preSamplerIterations)
        signals(sample-options.preSamplerIterations,:) = currentSignal;
    else
        presignals(sample,:) = currentSignal;
    end
    
%     if (mod(sample,10) == 1)
%         figure(1); clf; hold all;
%         plot(currentSignal);
%         plot(mean(signals(1:(sample-options.preSamplerIterations),:),1),'linewidth',2);
% % 
%         figure(2); clf;
%         plotParams(params);
%         q = [params(1:sample).muL];
%         r = reshape(cell2mat(q),numel(data.proxy),numel(q)/numel(data.proxy));
%         subplot(6,1,1); plot(r'); axis tight;
%         q = [params(1:sample).betaL];
%         r = reshape(cell2mat(q),numel(data.proxy),numel(q)/numel(data.proxy));
%         subplot(6,1,2); plot(r'); axis tight;
%         drawnow;        
%     end
%     pause;
%     fprintf('LNA:%d/%d\n',sample, totalIterations); 

    tracker('LNA',sample, totalIterations);

end

%% Store results

result.signals = signals;
result.presignals = presignals;
result.parameters = params;
result.times = data.target.times;
