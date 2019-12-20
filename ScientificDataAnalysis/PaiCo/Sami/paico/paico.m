function result = paico(data, options, tracker)
%PAICO Multiproxy climate reconstruction with pairwise comparisons
%   result = paico(data, options, tracker) calculates the climate signal
%   that best matches the pairwise comparisons of all proxy data. 
%
%   If options is omitted, defaults will be used. If tracker is omitted, no
%   output is generated.
%   
%   If data is set to 'info', result returns default options and
%   information about their meaning.
%
%   For information about data and result structures, type:
%   help rtstruct 


%   COPYRIGHT 
%       2012, Sami Hanhij�rvi <sami.hanhijarvi@iki.fi>
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
options = defaults(options);  

if(~exist('tracker','var'))
    % Empty function
    tracker = @(a,b,c)([]);
end

if (numel(options.regcov) == 1)
    options.regcov = eye(numel(data.target.times))/options.regcov;
else 
    options.regcov = inv(options.regcov);
end

if (options.resample.number > 0)
    if (~isfield(data.instrumental,'noisestd'))
        if (size(data.instrumental.data,1)>1)
            s = data.instrumental.data;
            s = s - repmat(nanmean(data.instrumental.data,1),size(data.instrumental.data,1),1);            
            data.instrumental.noisestd = sqrt(nanmean(s(:).^2));
        else
            error(['Instrumental noise STD was not set and could not be ' ...
                'calculated from given instrumental data. Set ' ...
                'data.instrumental.noisestd = 0 to ignore instrumental noise.']);
        end
    end
end




% Initialize comparison patterns
cache = inferProxies(data.target.times, data.proxy, tracker, options.fileCache);

% Select all pairwise comparisons information
[A counts] = selectProxies(cache, 1:cache.nrProxy, options.fileCache);    

% Find maximum likelihood solution
signal = optimize(A, counts, options);                

% Calibrate to instrumental data
result = calibrate(signal', data);

if (options.resample.number > 0)
    % Carry out resampling
    signals = zeros(options.resample.number, numel(signal));
    stds = zeros(options.resample.number, 1);
    if (options.resample.parallelize)        
        %ppm = ParforProgMon('Uncertainty ', options.resample.number, 1, 300, 80);
        %ppm = ParforProgMon('Uncertainty ', options.resample.number); CD: comment line
        parfor bi = 1:options.resample.number
            [A counts] = selectProxies(cache, randsample(cache.nrProxy, cache.nrProxy, true), options.fileCache);
            
            signal = optimize(A, counts, options);        

            br = calibrate(signal, data, true);
            signals(bi,:) = br.signal;
            stds(bi) = br.noisestd;
            %ppm.increment(); CD: comment line
        end
        %ppm.delete(); CD comment line
    else
        tracker('PAICO>resample');
        for bi = 1:options.resample.number
            [A counts] = selectProxies(cache, randsample(cache.nrProxy, cache.nrProxy, true), options.fileCache);
            
            signal = optimize(A, counts, options);        
            
            br = calibrate(signal, data, true);
            signals(bi,:) = br.signal;
            stds(bi) = br.noisestd;
            tracker('PAICO>resample', bi, options.resample.number);    
        end
    end
    
    result.resample.signals = signals;
    result.resample.noisestds = stds;
end