function result = barcast(data, options, tracker)
%BARCAST Multiproxy climate field reconstruction using Bayesian
%   Hierarchical models
%   result = barcast(data, options, tracker) calculates samples of climate 
%   field according to the model
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
%   Based on the code of Martin P. Tingley from Tingley et al. (2010).


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
    error('RECON:PAICO','At least one input parameter has to be defined.');    
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
options = defaults(options, data);  

if(~exist('tracker','var'))
    % Empty function
    tracker = @(a,b,c)([]);
end

%% Initialize

% Find which locations in proxies correspond to instrumental locations
for i = 1:numel(data.proxy)
    [~,ind] = min(earthDistances(data.instrumental.locations, data.proxy{i}.locations), [], 1);
    if (numel(ind)==1)
        ind = repmat(ind,1,size(data.proxy{i}.data,1));
    end
    data.proxy{i}.locationIndices = ind;
end
data.instrumental.locationIndices = 1:size(data.instrumental.locations,1);
% Add a single point to the timeline according to BARCAST. Assumes times are
% in years AD
timeline = data.target.times;


timeline = [min(timeline)-1 timeline];
model.timeline = sort(timeline,'ascend');

% Find which times in proxies correspond to timeline points
for i = 1:numel(data.proxy)
    [~, loc] = ismember(data.proxy{i}.times, model.timeline);
    data.proxy{i}.timeIndices = loc(loc>0);
end
[~, loc] = ismember(data.instrumental.times, model.timeline);
data.instrumental.timeIndices = loc(loc>0);

% Find distance matrix
model.distances = earthDistances(data.instrumental.locations);

% Draw initial values for the Bayesian model
[currentParams, currentField] = initialValues(data, model, options);

% Calculate spatial covariance
model.spatialCorrMatrix = exp(-currentParams.phi*model.distances);
model.invSpatialCorrMatrix = inv(model.spatialCorrMatrix);

if (~options.sampleCompleteField)
    % Discover missing patterns
    model.missingPatterns = findMissingPatterns(data, numel(model.timeline));

    % Calculate temporal covariance matrices for each missing pattern
    if (options.useSpatialCache)
        [model.spatialCovMatrices, model.sqrtSpatialCovMatrices] = calcSpatialCovariances(data, model, currentParams);
    end
end

params = currentParams;
params(options.samplerIterations).alpha = [];

tracker('BARCAST');
fields = zeros([size(currentField) options.samplerIterations]);

%% Sample MCMC chain
totalIterations = options.preSamplerIterations + options.samplerIterations;
for sample = 1:totalIterations;
    
    % Sample temperature field
    if (~options.sampleCompleteField)
        % Original Gibb's sampler from Tingley et al. slighly optimized.
        currentField(:,1) = sampleTemp0(model, currentParams, options.priors, currentField(:,2));
        for i = 2:size(currentField,2)-1
            currentField(:,i) = sampleTempk(data, model, currentParams, i, currentField(:,i-1), currentField(:,i+1));
        end
        currentField(:,end) = sampleTempLast(data, model, currentParams, currentField(:,end-1));
    else
        % Sample complete field at once
        currentField = sampleTempField(data, model, currentParams, options.priors);
    end
    % Sample autocorrelation coefficient
    currentParams.alpha = sampleAutocorrCoeff(model, currentParams, options.priors, currentField);
    
    % Sample autocorrelation mean parameter
    currentParams.mu = sampleAutocorrMean(model, currentParams, options.priors, currentField);
    
    % Sample spatial covariance spill parameter
    currentParams.sigma2 = sampleSpatialVariance(model, currentParams, options.priors, currentField);

    % Sample spatial covariance range parameter
    [currentParams.phi, model.spatialCorrMatrix, model.invSpatialCorrMatrix] = ...
        sampleSpatialCovarianceRange(model, currentParams, options.priors, options.MHpars, currentField);

    if (sample > options.preSamplerIterations)
        % Sample instrumental measurement error variance
        currentParams.tau2_I = sampleInstrumentalErrorVar(data, options.priors, currentField);
    end
    
    for i = 1:numel(data.proxy)
        % Sample proxy-specific parameters
        
        % Sample measurement error variance
        currentParams.tau2_P(i) = sampleProxyErrorVar(data, currentParams, options.priors, currentField, i);
        
        % Sample proxy multiplier parameter
        currentParams.Beta_1(i) = sampleProxyMultiplier(data, currentParams, options.priors, currentField, i);
        
        % Sample proxy multiplier parameter
        currentParams.Beta_0(i) = sampleProxyAddition(data, currentParams, options.priors, currentField, i);        
    end
            
    if (~options.sampleCompleteField)
        % Calculate spatial covariance matrices.
        if (options.useSpatialCache)    
            [model.spatialCovMatrices, model.sqrtSpatialCovMatrices] = ...
                calcSpatialCovariances(data, model, currentParams);    
        end
    end
    
    params(sample) = currentParams;                    
    if (sample > options.preSamplerIterations)

        fields(:,:,sample-options.preSamplerIterations) = currentField;
    end
%         
%     if (mod(sample,100) == 0)
%         figure(1); clf;
%         subplot(2,1,1);
%         imagesc(currentField(:,2:end));
%         subplot(2,1,2); hold all;
%         plot(mean(currentField(:,2:end),1));
%     %     
%         plot(data.target.times, mean(data.target.data,1));
%         plot(data.target.times, mean(data.proxy{1}.data,1));
% 
%         plot(data.instrumental.times, data.instrumental.data);
%         drawnow;
% 
%         figure(2);
%         plotParams(params);
%         drawnow;
%     end
% %     
%     pause;
%     fprintf('BARCAST:%d/%d\n',sample, totalIterations);


    tracker('BARCAST',sample, totalIterations);

end

%% Store results

result.fields = fields(:,2:end,:);
result.parameters = params;
result.times = timeline;
result.locations = data.instrumental.locations;
